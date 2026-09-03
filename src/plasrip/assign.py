#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Assign query plasmids to stable PTUs using the external PlasRiP database.

PlasRiP v0.2.x separates software and reference data. PTU assignment reads
``plasmid_reference_index.tsv`` and ``PTU_metadata.tsv`` from a versioned
reference database and uses the database's ``mash_index/`` sketches.
"""

from __future__ import annotations

import argparse
import os
import sys
import subprocess
from typing import List, Dict, Tuple, Optional

import pandas as pd
from tqdm import tqdm

from .database import (
    build_assignment_indexes,
    inspect_mash_sketch,
    load_database,
    ptu_metadata_lookup,
    validate_mash_runtime,
    validate_mash_sketch_compatibility,
    validate_query_mash_parameters,
    validate_reference_sketch_sample,
)
from .exceptions import PlasRipConfigError
from .schemas import PTU_COLUMNS, PTU_HITS_COLUMNS, reorder_columns

NETWORK_COVERAGE_CUTOFF = 0.6
MIN_NETWORK_MEMBERS = 3
# Keep representative scans portable without changing PTU-assignment semantics.
# The reference list is split only when the generated command would become large.
MASH_MULTI_MAX_COMMAND_CHARS = 24_000 if os.name == "nt" else 200_000
VERBOSE = True


def log(msg: str) -> None:
    if VERBOSE:
        sys.stderr.write(msg + "\n")


def ensure_sketch(
    fasta_path: str,
    mash_cmd: str,
    threads: int,
    mash_k: Optional[int] = None,
    mash_s: Optional[int] = None,
    mash_seed: Optional[int] = None,
    database=None,
) -> str:
    """Ensure a database-compatible Mash sketch exists for one query FASTA.

    Existing same-name ``.msh`` files are reused only after ``mash info``
    confirms that their hash function, seed, k-mer size, alphabet, sketch size,
    and number of sketches match the active reference database. Newly generated
    sketches explicitly set the database seed, k-mer size, and sketch size and
    are validated immediately after creation.
    """
    fasta_path = os.path.abspath(fasta_path)
    prefix, _ = os.path.splitext(fasta_path)
    msh_path = prefix + ".msh"

    if os.path.isfile(msh_path):
        if database is not None:
            info = inspect_mash_sketch(msh_path, mash_cmd)
            validate_mash_sketch_compatibility(
                database, info, label=f"Existing query sketch {msh_path}"
            )
        return msh_path

    cmd = [mash_cmd, "sketch"]
    if threads > 1:
        cmd += ["-p", str(threads)]
    if mash_seed is not None:
        cmd += ["-S", str(mash_seed)]
    if mash_k is not None:
        cmd += ["-k", str(mash_k)]
    if mash_s is not None:
        cmd += ["-s", str(mash_s)]
    cmd += ["-o", prefix, fasta_path]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except FileNotFoundError as e:
        raise PlasRipConfigError(f"Mash executable not found: {mash_cmd!r}") from e
    except subprocess.CalledProcessError as e:
        log(f"[ERROR] mash sketch failed: {' '.join(cmd)}")
        log(e.stderr)
        raise

    if database is not None:
        info = inspect_mash_sketch(msh_path, mash_cmd)
        validate_mash_sketch_compatibility(
            database, info, label=f"Generated query sketch {msh_path}"
        )
    return msh_path


def _reference_batches(
    ref_paths: List[str],
    max_command_chars: int = MASH_MULTI_MAX_COMMAND_CHARS,
) -> List[List[str]]:
    """Split reference paths into deterministic command-size batches.

    Mash accepts many query sketch files in one ``mash dist`` invocation, but
    passing thousands of absolute paths in a single command can exceed process
    command-line limits on some systems. Batching changes only process
    invocation size; result order and PTU-assignment logic are preserved.
    """
    if not ref_paths:
        return []

    batches: List[List[str]] = []
    batch: List[str] = []
    chars = 0
    for path in ref_paths:
        cost = len(str(path)) + 1
        if batch and chars + cost > max_command_chars:
            batches.append(batch)
            batch = []
            chars = 0
        batch.append(path)
        chars += cost
    if batch:
        batches.append(batch)
    return batches


def mash_dist_multi(
    query_msh: str,
    ref_paths: List[str],
    mash_cmd: str,
    threads: int
) -> List[str]:
    """Run ``mash dist`` for one query sketch against many reference sketches.

    Large reference lists are transparently processed in deterministic batches
    to avoid operating-system command-line length limits.
    """
    if not ref_paths:
        return []

    output_lines: List[str] = []
    batches = _reference_batches(ref_paths)
    for batch_no, batch in enumerate(batches, start=1):
        if threads > 1:
            cmd = [mash_cmd, "dist", "-p", str(threads), query_msh] + batch
        else:
            cmd = [mash_cmd, "dist", query_msh] + batch

        try:
            proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as e:
            head = " ".join(cmd[:5]) + (" ..." if len(cmd) > 5 else "")
            log(
                f"[ERROR] mash dist failed for representative batch "
                f"{batch_no}/{len(batches)}: {head}"
            )
            log(e.stderr)
            raise

        output_lines.extend(
            line for line in proc.stdout.splitlines() if line.strip()
        )

    return output_lines


def mash_dist_pair(
    query_msh: str,
    ref_msh: str,
    mash_cmd: str,
    threads: int
) -> Optional[float]:
    """
    Run `mash dist` of query vs a single ref, return distance or None.
    """
    if threads > 1:
        cmd = [mash_cmd, "dist", "-p", str(threads), query_msh, ref_msh]
    else:
        cmd = [mash_cmd, "dist", query_msh, ref_msh]

    try:
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        log(f"[ERROR] mash dist failed: {' '.join(cmd)}")
        log(e.stderr)
        return None

    for line in proc.stdout.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            return None
        try:
            return float(parts[2])
        except ValueError:
            return None

    return None


# ---------- parse representative hits ----------

def parse_rep_hits(
    fasta_path: str,
    lines: List[str],
    rep_key_to_record: Dict[str, Tuple[str, str]],
    threshold: float
) -> List[Dict[str, object]]:
    """Parse Mash hits against PTU representatives.

    ``rep_key_to_record`` maps Mash-internal reference IDs and stable reference
    aliases to ``(PTU_ID, representative_plasmid_id)``.
    """
    plasmid_name = os.path.basename(fasta_path)
    results: List[Dict[str, object]] = []

    for line in lines:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        ref_id = parts[1]
        try:
            dist = float(parts[2])
        except ValueError:
            continue
        if dist > threshold:
            continue

        ref_base = os.path.basename(ref_id)
        ref_stem, _ = os.path.splitext(ref_base)
        rec = (
            rep_key_to_record.get(ref_id)
            or rep_key_to_record.get(ref_base)
            or rep_key_to_record.get(ref_stem)
        )
        if rec is None:
            continue
        ptu_id, representative_id = rec
        results.append({
            "plasmid": plasmid_name,
            "PTU_ID": ptu_id,
            "Representative": representative_id,
            "dist": dist,
        })

    return results


# ---------- full-member + network stats for a PTU ----------

def full_member_check(
    query_msh: str,
    member_paths: List[str],
    threshold: float,
    mash_cmd: str,
    threads: int
) -> Tuple[bool, Optional[float], Optional[float], int,
           Optional[float], int, Optional[float]]:
    """
    For one PTU: compare query vs ALL member .msh files.

    Returns:
      (all_pass, max_dist_all, avg_dist_all, n_members,
       coverage, n_within, mean_dist_within)
    """
    if not member_paths:
        return False, None, None, 0, None, 0, None

    dists: List[float] = []
    for mp in member_paths:
        d = mash_dist_pair(query_msh, mp, mash_cmd, threads)
        if d is None:
            # treat as failure: we do not trust this PTU for full-membership
            return False, None, None, len(dists), None, 0, None
        dists.append(d)

    if not dists:
        return False, None, None, 0, None, 0, None

    n_members = len(dists)
    max_d = max(dists)
    avg_d = sum(dists) / n_members

    within = [d for d in dists if d <= threshold]
    n_within = len(within)
    coverage = n_within / n_members if n_members > 0 else None
    mean_in = sum(within) / n_within if n_within > 0 else None

    all_pass = (n_within == n_members and n_members > 0)

    return all_pass, max_d, avg_d, n_members, coverage, n_within, mean_in


# ---------- PTU assignment ----------

def assign_ptu(
    plasmid_name: str,
    query_msh: Optional[str],
    hits: List[Dict[str, object]],
    ptu_to_members: Dict[str, List[str]],
    threshold: float,
    mash_cmd: str,
    threads: int,
) -> Dict[str, object]:
    """
    Assign PTU for one plasmid; return a summary dict.
    """
    # default: no_hit
    res = {
        "query": plasmid_name,
        "n_ptu_hits": 0,
        "ptu_assigned": "NA",
        "how_assigned": "no_hit",
        "closest_ptu": "NA",
        "closest_rep": "NA",
        "closest_rep_dist": "NA",
        "second_ptu": "NA",
        "second_rep_dist": "NA",
        "rep_dist_gap": "NA",
        "all_member_ptu": "NA",
        "all_member_max_dist": "NA",
        "all_member_avg_dist": "NA",
        "all_member_n": "NA",
        "network_ptu": "NA",
        "network_coverage": "NA",
        "network_mean_dist_leq_threshold": "NA",
        "network_n_within_threshold": "NA",
    }

    if query_msh is None:
        res["how_assigned"] = "failed"
        return res

    if not hits:
        return res

    # PTU -> best representative distance
    ptu_info: Dict[str, Dict[str, object]] = {}
    for h in hits:
        ptu_id = str(h["PTU_ID"])
        d = float(h["dist"])
        rep = str(h["Representative"])
        info = ptu_info.get(ptu_id)
        if info is None or d < info["best_rep_dist"]:
            ptu_info[ptu_id] = {
                "best_rep": rep,
                "best_rep_dist": d,
            }

    n_candidates = len(ptu_info)
    res["n_ptu_hits"] = n_candidates

    # sort PTUs by best representative distance
    sorted_ptus = sorted(ptu_info.items(), key=lambda x: x[1]["best_rep_dist"])
    best_ptu, best_info = sorted_ptus[0]
    best_rep = best_info["best_rep"]
    best_dist = best_info["best_rep_dist"]

    if len(sorted_ptus) > 1:
        second_ptu, second_info = sorted_ptus[1]
        second_dist = second_info["best_rep_dist"]
        gap = second_dist - best_dist
    else:
        second_ptu, second_dist, gap = None, None, None

    res["closest_ptu"] = best_ptu
    res["closest_rep"] = best_rep
    res["closest_rep_dist"] = best_dist
    if second_ptu is not None:
        res["second_ptu"] = second_ptu
        res["second_rep_dist"] = second_dist
        res["rep_dist_gap"] = gap

    # ---- compute member-based stats for all candidate PTUs ----
    member_stats: Dict[str, Dict[str, object]] = {}
    for ptu_id in ptu_info.keys():
        member_paths = ptu_to_members.get(ptu_id, [])
        if not member_paths:
            continue
        (all_pass, max_d, avg_d, n_m,
         cov, n_within, mean_in) = full_member_check(
            query_msh, member_paths, threshold, mash_cmd, threads
        )
        member_stats[ptu_id] = {
            "all_pass": all_pass,
            "max_dist": max_d,
            "avg_dist": avg_d,
            "n_members": n_m,
            "coverage": cov,
            "n_within": n_within,
            "mean_in": mean_in,
        }

    # ---- full-member PTUs ----
    passing = [
        (ptu_id, st)
        for ptu_id, st in member_stats.items()
        if st["all_pass"]
    ]

    if passing:
        # among full-pass PTUs: min(max_dist, then avg_dist, then rep_dist)
        def key_fn(item):
            ptu_id, stat = item
            max_d = stat["max_dist"] if stat["max_dist"] is not None else 1e9
            avg_d = stat["avg_dist"] if stat["avg_dist"] is not None else 1e9
            rep_d = ptu_info[ptu_id]["best_rep_dist"]
            return (max_d, avg_d, rep_d)

        passing_sorted = sorted(passing, key=key_fn)
        chosen_ptu, stat = passing_sorted[0]

        res["ptu_assigned"] = chosen_ptu
        res["how_assigned"] = (
            "all_members" if len(passing) == 1 else "all_members_multi"
        )
        res["all_member_ptu"] = chosen_ptu
        res["all_member_max_dist"] = stat["max_dist"]
        res["all_member_avg_dist"] = stat["avg_dist"]
        res["all_member_n"] = stat["n_members"]

        # network summary for the chosen full-member PTU
        res["network_ptu"] = chosen_ptu
        res["network_coverage"] = stat["coverage"]
        res["network_mean_dist_leq_threshold"] = stat["mean_in"]
        res["network_n_within_threshold"] = stat["n_within"]

        return res

    # ---- no full-member PTU; compute best network support ----
    network_candidates = [
        (ptu_id, st)
        for ptu_id, st in member_stats.items()
        if st["coverage"] is not None and st["n_members"] > 0 and st["n_within"] > 0
    ]

    if network_candidates:
        def net_key(item):
            ptu_id, stat = item
            cov = stat["coverage"]
            # tie-breaker: smaller representative distance
            rep_d = ptu_info[ptu_id]["best_rep_dist"]
            return (cov, -rep_d)

        best_net_ptu, best_net_stat = max(network_candidates, key=net_key)
        res["network_ptu"] = best_net_ptu
        res["network_coverage"] = best_net_stat["coverage"]
        res["network_mean_dist_leq_threshold"] = best_net_stat["mean_in"]
        res["network_n_within_threshold"] = best_net_stat["n_within"]
    else:
        best_net_ptu = None
        best_net_stat = None

    # ---- if only one candidate PTU: keep behaviour simple ----
    if n_candidates == 1:
        res["ptu_assigned"] = best_ptu
        res["how_assigned"] = "single_ptu_by_rep"
        return res

    # ---- multiple PTUs, no full-member; try network_support ----
    if (
        best_net_ptu is not None
        and best_net_stat is not None
        and best_net_stat["coverage"] is not None
        and best_net_stat["n_members"] >= MIN_NETWORK_MEMBERS
    ):
        cov = best_net_stat["coverage"]
        rep_d_net = ptu_info[best_net_ptu]["best_rep_dist"]
        if cov >= NETWORK_COVERAGE_CUTOFF and rep_d_net <= threshold:
            res["ptu_assigned"] = best_net_ptu
            res["how_assigned"] = "network_support"
            return res

    # ---- fallback: closest representative + distance gap logic ----
    gap_cutoff = 0.01

    if best_dist <= threshold:
        if second_ptu is None or gap is None or gap >= gap_cutoff:
            res["ptu_assigned"] = best_ptu
            res["how_assigned"] = "closest_rep"
        else:
            res["ptu_assigned"] = "NA"
            res["how_assigned"] = "ambiguous"
    else:
        res["ptu_assigned"] = "NA"
        res["how_assigned"] = "above_threshold"

    return res


# ---------- public output formatting ----------

NO_PTU_VALUE = "NA"
NO_PTU_PGG_STATUS = "NOT_APPLICABLE_NO_PTU"


def _normalise_query_id(value: object) -> str:
    """Return the stable query identifier used in public TSV outputs."""
    base = os.path.basename(str(value))
    stem, _ = os.path.splitext(base)
    return stem


def _format_assignment_output(assignments: List[Dict[str, object]]) -> pd.DataFrame:
    """Convert internal assignment records to the frozen v0.2.0 public schema.

    PGG semantics:
      * assigned PTU + stable PGG: report PGG ID and database assignment status;
      * assigned PTU without stable PGG: leave PGG ID empty and retain the exact
        database unassigned status;
      * no PTU assignment: report ``PGG_ID=NA`` and
        ``PGG_assignment_status=NOT_APPLICABLE_NO_PTU`` because no PTU exists
        from which a PGG assignment could be defined.
    """
    rows: List[Dict[str, object]] = []
    for rec in assignments:
        ptu = str(rec.get("ptu_assigned", NO_PTU_VALUE) or NO_PTU_VALUE)
        if ptu == NO_PTU_VALUE:
            pgg_id = NO_PTU_VALUE
            pgg_status = NO_PTU_PGG_STATUS
        else:
            pgg_id = str(rec.get("PGG_ID", "") or "")
            pgg_status = str(rec.get("PGG_assignment_status", "") or "")

        rows.append({
            "PLASMID": _normalise_query_id(rec.get("query", "")),
            "PTU": ptu,
            "PGG_ID": pgg_id,
            "PGG_assignment_status": pgg_status,
            "PTU_hits": rec.get("n_ptu_hits", 0),
            "how_assigned": rec.get("how_assigned", ""),
            "closest_rep_dist": rec.get("closest_rep_dist", NO_PTU_VALUE),
            "PTU_closest": rec.get("closest_ptu", NO_PTU_VALUE),
            "PTU_closest_rep": rec.get("closest_rep", NO_PTU_VALUE),
            "PTU_second": rec.get("second_ptu", NO_PTU_VALUE),
            "PTU_second_rep_dist": rec.get("second_rep_dist", NO_PTU_VALUE),
            "PTU_rep_dist_gap": rec.get("rep_dist_gap", NO_PTU_VALUE),
            "PTU_all_member": rec.get("all_member_ptu", NO_PTU_VALUE),
            "PTU_all_member_max_dist": rec.get("all_member_max_dist", NO_PTU_VALUE),
            "PTU_all_member_avg_dist": rec.get("all_member_avg_dist", NO_PTU_VALUE),
            "PTU_all_member_n": rec.get("all_member_n", NO_PTU_VALUE),
            "PTU_network": rec.get("network_ptu", NO_PTU_VALUE),
            "PTU_network_coverage": rec.get("network_coverage", NO_PTU_VALUE),
            "PTU_network_mean_dist_leq_threshold": rec.get(
                "network_mean_dist_leq_threshold", NO_PTU_VALUE
            ),
            "PTU_network_n_within_threshold": rec.get(
                "network_n_within_threshold", NO_PTU_VALUE
            ),
        })

    return reorder_columns(pd.DataFrame(rows), PTU_COLUMNS, keep_extra=False)


def _format_hits_output(hits: List[Dict[str, object]]) -> pd.DataFrame:
    """Convert representative-hit records to the frozen public hit schema."""
    rows: List[Dict[str, object]] = []
    for hit in hits:
        rows.append({
            "PLASMID": _normalise_query_id(hit.get("plasmid", "")),
            "PTU_ID": hit.get("PTU_ID", ""),
            "PGG_ID": hit.get("PGG_ID", ""),
            "PGG_assignment_status": hit.get("PGG_assignment_status", ""),
            "Representative": hit.get("Representative", ""),
            "dist": hit.get("dist", ""),
        })
    return reorder_columns(pd.DataFrame(rows), PTU_HITS_COLUMNS, keep_extra=False)


# ---------- high-level runner ----------

def run_ptu_assignment(
    input_fastas: List[str],
    database_dir: str | None,
    output: str,
    threshold: float | None = None,
    threads: int = 1,
    mash_cmd: str = "mash",
    mash_k: Optional[int] = None,
    mash_s: Optional[int] = None,
    show_progress: bool = True,
    _database=None,
    _assignment_indexes=None,
    _ptu_meta_lookup=None,
) -> pd.DataFrame:
    """Assign query plasmids to PTUs using a PlasRiP reference database."""
    global VERBOSE
    VERBOSE = show_progress

    db = _database or load_database(database_dir, check_sketches=True)
    runtime_mash_version = validate_mash_runtime(db, mash_cmd)
    checked_reference_sketches = validate_reference_sketch_sample(
        db, mash_cmd, sample_size=3
    )
    db_k, db_s = validate_query_mash_parameters(db, mash_k, mash_s)
    cutoff = db.ptu_threshold if threshold is None else float(threshold)
    if not 0.0 <= cutoff <= 1.0:
        raise PlasRipConfigError(
            f"PTU Mash distance cutoff must be within [0, 1], got {cutoff}."
        )
    if int(threads) < 1:
        raise PlasRipConfigError(f"threads must be >= 1, got {threads}.")

    log("========== PlasRiP PTU assignment ==========")
    log(f"[DB  ] root          : {db.root}")
    log(f"[DB  ] version       : {db.database_version}")
    log(f"[DB  ] mash index    : {db.mash_index_dir}")
    log(f"[RUN ] mash          : {mash_cmd} (version {runtime_mash_version})")
    log(f"[RUN ] threshold     : {cutoff}")
    log(f"[RUN ] mash -S/-k/-s : {db.mash_seed}/{db_k}/{db_s}")
    log(f"[RUN ] reference QC  : {len(checked_reference_sketches)} representative sketches")
    log(f"[RUN ] mash threads  : {threads}")

    input_paths: List[str] = []
    for p in input_fastas:
        if not os.path.isfile(p):
            raise FileNotFoundError(p)
        input_paths.append(os.path.abspath(p))

    if _assignment_indexes is None:
        ptu_to_members, rep_paths, rep_key_to_record = build_assignment_indexes(db)
    else:
        ptu_to_members, rep_paths, rep_key_to_record = _assignment_indexes
    ptu_meta = _ptu_meta_lookup or ptu_metadata_lookup(db)
    log(f"[RUN ] plasmids       : {len(input_paths)} files")
    log(f"[RUN ] representatives : {len(rep_paths)} sketches")

    all_hits: List[Dict[str, object]] = []
    hits_by_plasmid: Dict[str, List[Dict[str, object]]] = {}
    query_msh_map: Dict[str, Optional[str]] = {}

    iterator = (
        tqdm(input_paths, total=len(input_paths), desc="Plasmids", unit="file", ncols=100)
        if show_progress else input_paths
    )

    for fa in iterator:
        base = os.path.basename(fa)
        if show_progress and hasattr(iterator, "set_postfix_str"):
            iterator.set_postfix_str(base)
        try:
            query_msh = ensure_sketch(
                fa, mash_cmd, threads, mash_k=db_k, mash_s=db_s,
                mash_seed=db.mash_seed, database=db
            )
        except PlasRipConfigError:
            # Parameter/version incompatibilities affect scientific validity and
            # must abort the run rather than silently marking one query failed.
            raise
        except Exception:
            log(f"[ERROR] sketch failed, skip: {base}")
            query_msh_map[base] = None
            hits_by_plasmid[base] = []
            continue

        query_msh_map[base] = query_msh
        try:
            lines = mash_dist_multi(query_msh, rep_paths, mash_cmd, threads)
        except Exception:
            log(f"[ERROR] mash dist (representatives) failed, skip: {base}")
            hits_by_plasmid[base] = []
            continue

        hits = parse_rep_hits(fa, lines, rep_key_to_record, cutoff)
        for hit in hits:
            meta = ptu_meta.get(str(hit["PTU_ID"]), {})
            hit["PGG_ID"] = meta.get("PGG_ID", "")
            hit["PGG_assignment_status"] = meta.get("PGG_assignment_status", "")
        hits_by_plasmid[base] = hits
        all_hits.extend(hits)
        log(f"[HITS] {base}: {len(hits)} PTU(s) with dist <= {cutoff}")

    assignments: List[Dict[str, object]] = []
    for fa in input_paths:
        base = os.path.basename(fa)
        res = assign_ptu(
            plasmid_name=base,
            query_msh=query_msh_map.get(base),
            hits=hits_by_plasmid.get(base, []),
            ptu_to_members=ptu_to_members,
            threshold=cutoff,
            mash_cmd=mash_cmd,
            threads=threads,
        )
        ptu_id = str(res.get("ptu_assigned", "NA"))
        if ptu_id and ptu_id != "NA":
            meta = ptu_meta.get(ptu_id, {})
            res["PGG_ID"] = meta.get("PGG_ID", "")
            res["PGG_assignment_status"] = meta.get("PGG_assignment_status", "")
        else:
            res["PGG_ID"] = ""
            res["PGG_assignment_status"] = ""
        assignments.append(res)
        log(
            f"[ASSN] {base}: PTU={res['ptu_assigned']}, "
            f"PGG={res.get('PGG_ID') or 'NA'}, method={res['how_assigned']}, "
            f"hits={res['n_ptu_hits']}, closest_dist={res['closest_rep_dist']}"
        )

    df_assign = _format_assignment_output(assignments)
    df_hits = _format_hits_output(all_hits)

    hits_path = os.path.abspath(output) + ".hits.tsv"
    # Always write the candidate-hit table, including a header-only file when
    # no representative PTU passes the Mash cutoff. This gives downstream
    # workflows a stable two-file output contract.
    df_hits.to_csv(hits_path, sep="\t", index=False)

    out_path = os.path.abspath(output)
    df_assign.to_csv(out_path, sep="\t", index=False)

    n_total = len(df_assign)
    n_assigned = (df_assign["PTU"] != NO_PTU_VALUE).sum() if n_total else 0
    n_pgg = (
        ((df_assign["PGG_ID"] != NO_PTU_VALUE) & (df_assign["PGG_ID"] != "")).sum()
        if n_total else 0
    )
    log("=============================================")
    log(f"[DONE] summary table    : {out_path}")
    log(f"[DONE] candidate hits   : {hits_path}")
    log(f"[DONE] total plasmids   : {n_total}")
    log(f"[DONE] assigned to PTU  : {n_assigned}")
    log(f"[DONE] assigned PGG     : {n_pgg}")
    return df_assign


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Assign plasmids to PTUs using the PlasRiP reference database"
    )
    p.add_argument("-i", "--input", nargs="+", required=True, help="Input plasmid FASTA files")
    p.add_argument(
        "--database",
        default=None,
        help="Reference database directory (or set PLASRIP_DATABASE)",
    )
    p.add_argument("-o", "--output", required=True, help="Output PTU summary TSV")
    p.add_argument(
        "--ptu-threshold",
        type=float,
        default=None,
        help="PTU Mash cutoff; default is read from database_manifest.tsv",
    )
    p.add_argument("-t", "--threads", type=int, default=1)
    p.add_argument("--mash", default="mash")
    p.add_argument(
        "--mash-k",
        type=int,
        default=None,
        help="Advanced: must match the reference database manifest when supplied",
    )
    p.add_argument(
        "--mash-sketch-size",
        type=int,
        default=None,
        help="Advanced: must match the reference database manifest when supplied",
    )
    p.add_argument("--quiet", action="store_true")
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)
    run_ptu_assignment(
        input_fastas=args.input,
        database_dir=args.database,
        output=args.output,
        threshold=args.ptu_threshold,
        threads=args.threads,
        mash_cmd=args.mash,
        mash_k=args.mash_k,
        mash_s=args.mash_sketch_size,
        show_progress=not args.quiet,
    )


if __name__ == "__main__":
    main()

run_assign = run_ptu_assignment
