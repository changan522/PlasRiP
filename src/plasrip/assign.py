#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PlasRip module: assign plasmids to PTUs using Mash distance
and an existing PTU member library.

Public API:
    run_ptu_assignment(
        input_fastas, rep_dir, output,
        threshold=0.05,
        threads=1,
        mash_cmd="mash",
        mash_k=None,
        mash_s=None,
        show_progress=True,
    )
    main()  # CLI entry

Also exported for reuse in the main pipeline (to avoid re-reading big tables):
    read_csv_from_data(...)
    build_member_index(...)
    build_rep_index(...)
    ensure_sketch(...)
    mash_dist_multi(...)
    parse_rep_hits(...)
    assign_ptu(...)

Inputs (fixed file names, inside package data directory):
  - PTU-metadata.csv  (columns: PTU_ID, Representative [.msh name])
  - PTUmember.csv     (columns: Plasmid [.msh name], PTU)
"""

import argparse
import os
import sys
import subprocess
from typing import List, Dict, Tuple, Optional

import pandas as pd
from tqdm import tqdm

# Fixed file names (relative to this module's data/)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "data")
PTU_METADATA_FILENAME = "PTU-metadata.csv"
PTU_MEMBER_FILENAME = "PTUmember.csv"

# network support cutoff (for multi-PTU hits, when no full-membership PTU)
NETWORK_COVERAGE_CUTOFF = 0.6   # require at least 60% of members within threshold
MIN_NETWORK_MEMBERS = 3         # only trust network support if PTU has >= this many members


# ---------- basic utils ----------

VERBOSE = True

def log(msg: str) -> None:
    """Unified log to stderr when verbose mode is enabled."""
    if VERBOSE:
        sys.stderr.write(msg + "\n")


def read_csv_from_data(
    filename: str,
    required_cols: List[str],
    label: str
) -> pd.DataFrame:
    """
    Read a CSV/TSV file from plasrip/data, ensure required columns exist.
    Auto-detect separator, handle BOM.
    """
    path = os.path.join(DATA_DIR, filename)
    if not os.path.isfile(path):
        log(f"[ERROR] {filename} not found in data directory.")
        log(f"        Please place {filename} under: {DATA_DIR}")
        raise FileNotFoundError(path)

    df = pd.read_csv(path, encoding="utf-8-sig", sep=None, engine="python")
    df.columns = df.columns.str.replace("\ufeff", "").str.strip()

    missing = set(required_cols) - set(df.columns)
    if missing:
        raise ValueError(
            f"[ERROR] {label} file {path} is missing columns: {missing}. "
            f"Found columns: {list(df.columns)}"
        )

    log(f"[SETUP] Loaded {label} file: {path} "
        f"({len(df)} rows)")
    return df[required_cols].copy()


# ---------- build PTU index ----------

def build_member_index(df_member: pd.DataFrame, member_dir: str) -> Dict[str, List[str]]:
    """
    Build:
      ptu_to_members[PTU_ID] = [absolute paths of member .msh files]
    """
    member_dir = os.path.abspath(member_dir)
    ptu_to_members: Dict[str, List[str]] = {}
    missing: List[str] = []

    for _, row in df_member.iterrows():
        ptu_id = str(row["PTU"])
        plasmid_name = os.path.basename(str(row["Plasmid"]).strip())
        path = os.path.join(member_dir, plasmid_name)

        if not os.path.isfile(path):
            missing.append(path)
            continue

        ptu_to_members.setdefault(ptu_id, []).append(path)

    log(
        f"[SETUP] Member index: {len(df_member)} member rows, "
        f"{len(ptu_to_members)} PTUs with .msh members"
    )
    if missing:
        log(f"[WARN ] Missing member .msh files: {len(missing)} (show up to 10)")
        for p in missing[:10]:
            log(f"        - {p}")

    return ptu_to_members


def build_rep_index(
    df_meta: pd.DataFrame,
    rep_dir: str
) -> Tuple[List[str], Dict[str, str]]:
    """
    Build representative library from PTU-metadata:
      - rep_paths: list of representative .msh paths
      - rep_key_to_ptu: ref-name keys -> PTU_ID

    For Representative = GCA_xxx_genomic_1.msh we create keys:
      - "GCA_xxx_genomic_1.msh"
      - "GCA_xxx_genomic_1"
    """
    rep_dir = os.path.abspath(rep_dir)
    rep_paths: List[str] = []
    rep_key_to_ptu: Dict[str, str] = {}
    missing: List[str] = []

    for _, row in df_meta.iterrows():
        ptu_id = str(row["PTU_ID"])
        rep_name = os.path.basename(str(row["Representative"]).strip())
        stem, _ = os.path.splitext(rep_name)
        rep_path = os.path.join(rep_dir, rep_name)

        if os.path.isfile(rep_path):
            rep_paths.append(rep_path)
        else:
            missing.append(rep_path)

        # multiple string forms mapping to the same PTU
        for key in (rep_name, stem):
            if key in rep_key_to_ptu and rep_key_to_ptu[key] != ptu_id:
                log(
                    f"[WARN ] Representative key {key} "
                    f"mapped to PTU {rep_key_to_ptu[key]} and {ptu_id}; keep first"
                )
                continue
            rep_key_to_ptu[key] = ptu_id

    log(
        f"[SETUP] Rep index: {len(df_meta)} PTUs, "
        f"{len(rep_paths)} representative .msh, "
        f"{len(rep_key_to_ptu)} name keys"
    )
    if missing:
        log(f"[WARN ] Missing representative .msh: {len(missing)} (show up to 10)")
        for p in missing[:10]:
            log(f"        - {p}")

    return rep_paths, rep_key_to_ptu


# ---------- Mash wrappers ----------

def ensure_sketch(
    fasta_path: str,
    mash_cmd: str,
    threads: int,
    mash_k: Optional[int] = None,
    mash_s: Optional[int] = None,
) -> str:
    """
    Ensure a .msh exists for given FASTA:
      - if same-name .msh exists: reuse
      - else: run `mash sketch` (with -p threads if >1, plus optional -k/-s)
    """
    fasta_path = os.path.abspath(fasta_path)
    prefix, _ = os.path.splitext(fasta_path)
    msh_path = prefix + ".msh"

    if os.path.isfile(msh_path):
        return msh_path

    cmd = [mash_cmd, "sketch"]
    if threads > 1:
        cmd += ["-p", str(threads)]
    if mash_k is not None:
        cmd += ["-k", str(mash_k)]
    if mash_s is not None:
        cmd += ["-s", str(mash_s)]
    cmd += [fasta_path, "-o", prefix]

    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        log(f"[ERROR] mash sketch failed: {' '.join(cmd)}")
        log(e.stderr)
        raise

    return msh_path


def mash_dist_multi(
    query_msh: str,
    ref_paths: List[str],
    mash_cmd: str,
    threads: int
) -> List[str]:
    """
    Run `mash dist` of query vs many refs, return stdout lines.
    """
    if threads > 1:
        cmd = [mash_cmd, "dist", "-p", str(threads), query_msh] + ref_paths
    else:
        cmd = [mash_cmd, "dist", query_msh] + ref_paths

    try:
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        head = " ".join(cmd[:5]) + (" ..." if len(cmd) > 5 else "")
        log(f"[ERROR] mash dist failed: {head}")
        log(e.stderr)
        raise

    return [l for l in proc.stdout.strip().split("\n") if l.strip()]


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
    rep_key_to_ptu: Dict[str, str],
    threshold: float
) -> List[Dict[str, object]]:
    """
    Parse mash dist output vs representatives for a single query.

    Keep only dist <= threshold and successfully mapped PTU.

    Returns records:
      plasmid, PTU_ID, Representative(stem), dist
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

        ptu_id = (
            rep_key_to_ptu.get(ref_stem)
            or rep_key_to_ptu.get(ref_base)
            or rep_key_to_ptu.get(ref_id)
        )
        if ptu_id is None:
            continue

        results.append({
            "plasmid": plasmid_name,
            "PTU_ID": ptu_id,
            "Representative": ref_stem,
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


# ---------- high-level runner ----------

def run_ptu_assignment(
    input_fastas: List[str],
    rep_dir: str,
    output: str,
    threshold: float = 0.05,
    threads: int = 1,
    mash_cmd: str = "mash",
    mash_k: Optional[int] = None,
    mash_s: Optional[int] = None,
    show_progress: bool = True,
) -> pd.DataFrame:
    """
    High-level function: run PTU assignment for multiple plasmid FASTA files.
    """
    global VERBOSE
    VERBOSE = show_progress
    log("========== PlasRip assignment ==========")
    log(f"[RUN ] mash          : {mash_cmd}")
    log(f"[RUN ] threshold     : {threshold}")
    log(f"[RUN ] mash threads  : {threads}")
    if mash_k is not None:
        log(f"[RUN ] mash -k       : {mash_k}")
    if mash_s is not None:
        log(f"[RUN ] mash -s       : {mash_s}")

    # check input FASTA files
    input_paths: List[str] = []
    for p in input_fastas:
        if not os.path.isfile(p):
            log(f"[ERROR] Input FASTA not found: {p}")
            raise FileNotFoundError(p)
        input_paths.append(os.path.abspath(p))
    log(f"[RUN ] plasmids       : {len(input_paths)} files")

    # load PTU metadata & members from package data
    df_meta = read_csv_from_data(
        PTU_METADATA_FILENAME,
        ["PTU_ID", "Representative"],
        label="PTU metadata"
    )
    df_member = read_csv_from_data(
        PTU_MEMBER_FILENAME,
        ["Plasmid", "PTU"],
        label="PTU members"
    )

    # build indexes
    ptu_to_members = build_member_index(df_member, rep_dir)
    rep_paths, rep_key_to_ptu = build_rep_index(df_meta, rep_dir)
    log(f"[RUN ] representatives : {len(rep_paths)} .msh files")
    log("========================================")

    all_hits: List[Dict[str, object]] = []
    hits_by_plasmid: Dict[str, List[Dict[str, object]]] = {}
    query_msh_map: Dict[str, Optional[str]] = {}

    # Step 1: sketch + representative scan
    if show_progress:
        iterator = tqdm(
            input_paths,
            total=len(input_paths),
            desc="Plasmids",
            unit="file",
            ncols=100,
        )
    else:
        iterator = input_paths

    for fa in iterator:
        base = os.path.basename(fa)
        if show_progress and hasattr(iterator, "set_postfix_str"):
            iterator.set_postfix_str(base)

        # sketch
        try:
            query_msh = ensure_sketch(
                fa, mash_cmd, threads, mash_k=mash_k, mash_s=mash_s
            )
        except Exception:
            log(f"[ERROR] sketch failed, skip: {base}")
            query_msh_map[base] = None
            hits_by_plasmid[base] = []
            continue

        query_msh_map[base] = query_msh

        # representative mash dist
        try:
            lines = mash_dist_multi(query_msh, rep_paths, mash_cmd, threads)
        except Exception:
            log(f"[ERROR] mash dist (reps) failed, skip: {base}")
            hits_by_plasmid[base] = []
            continue

        hits = parse_rep_hits(fa, lines, rep_key_to_ptu, threshold)
        hits_by_plasmid[base] = hits
        all_hits.extend(hits)

        log(f"[HITS] {base}: {len(hits)} PTU(s) with dist <= {threshold}")

    # Step 2: assign PTU per plasmid
    assignments: List[Dict[str, object]] = []

    for fa in input_paths:
        base = os.path.basename(fa)
        query_msh = query_msh_map.get(base)
        hits = hits_by_plasmid.get(base, [])
        res = assign_ptu(
            plasmid_name=base,
            query_msh=query_msh,
            hits=hits,
            ptu_to_members=ptu_to_members,
            threshold=threshold,
            mash_cmd=mash_cmd,
            threads=threads,
        )
        assignments.append(res)
        log(
            f"[ASSN] {base}: PTU={res['ptu_assigned']}, "
            f"method={res['how_assigned']}, "
            f"hits={res['n_ptu_hits']}, "
            f"closest_dist={res['closest_rep_dist']}, "
            f"gap={res['rep_dist_gap']}"
        )

    df_assign = pd.DataFrame(assignments)

    # write hits table
    hits_path = os.path.abspath(output) + ".hits.tsv"
    if all_hits:
        pd.DataFrame(all_hits).to_csv(hits_path, sep="\t", index=False)
        log(f"[OUT ] rep hits table  : {hits_path}")
    else:
        log("[OUT ] rep hits table  : empty (no hits)")

    # write summary table
    out_path = os.path.abspath(output)
    df_assign.to_csv(out_path, sep="\t", index=False)

    # small summary
    n_total = len(df_assign)
    n_assigned = (df_assign["ptu_assigned"] != "NA").sum()
    n_amb = (df_assign["how_assigned"] == "ambiguous").sum()
    n_nohit = (df_assign["how_assigned"] == "no_hit").sum()

    log("========================================")
    log(f"[DONE] summary table    : {out_path}")
    log(f"[DONE] total plasmids   : {n_total}")
    log(f"[DONE] assigned to PTU  : {n_assigned}")
    log(f"[DONE] ambiguous        : {n_amb}")
    log(f"[DONE] no_hit           : {n_nohit}")
    log("========================================")

    return df_assign


# ---------- CLI ----------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Assign plasmids to PTUs using Mash and PTU member library"
    )
    p.add_argument(
        "-i", "--input",
        nargs="+",
        required=True,
        help="Input plasmid FASTA files, e.g. -i p1.fasta p2.fasta"
    )
    p.add_argument(
        "-r", "--rep-dir",
        required=True,
        help="Directory containing all PTU member .msh files (including reps)"
    )
    p.add_argument(
        "-o", "--output",
        required=True,
        help="Output summary TSV (one row per plasmid)"
    )
    p.add_argument(
        "--threshold",
        type=float,
        default=0.05,
        help="Mash distance cutoff (default: 0.05)"
    )
    p.add_argument(
        "-t", "--threads",
        type=int,
        default=1,
        help="Threads for mash (-p), default: 1"
    )
    p.add_argument(
        "--mash",
        default="mash",
        help="Mash executable (default: mash)"
    )
    p.add_argument(
        "--mash-k",
        type=int,
        default=None,
        help="Mash sketch -k (k-mer size). Default: Mash internal default (21)."
    )
    p.add_argument(
        "--mash-sketch-size",
        type=int,
        default=None,
        help="Mash sketch -s (sketch size). Default: Mash internal default (1000)."
    )
    return p.parse_args()


def main() -> None:
    """CLI entry for `PlasRip-assign`."""
    args = parse_args()
    run_ptu_assignment(
        input_fastas=args.input,
        rep_dir=args.rep_dir,
        output=args.output,
        threshold=args.threshold,
        threads=args.threads,
        mash_cmd=args.mash,
        mash_k=args.mash_k,
        mash_s=args.mash_sketch_size,
        show_progress=True,
    )


if __name__ == "__main__":
    main()
    
run_assign = run_ptu_assignment