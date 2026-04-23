#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PlasRip submodule: plasmid mobility prediction."""

from __future__ import annotations

import os
import subprocess
import sys
from typing import Any, Dict, Iterable, List, Optional

import pandas as pd

def _msg(msg: str, show_progress: bool = True, is_err: bool = False) -> None:
    if not show_progress:
        return
    stream = sys.stderr
    stream.write(msg + "\n")
    stream.flush()

def _run_cmd(cmd: List[str], show_progress: bool = True) -> None:
    _msg(f"[mobility] RUN: {' '.join(cmd)}", show_progress=show_progress)
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError as e:
        _msg(f"[mobility] ERROR: command failed (exit code={e.returncode})", show_progress=show_progress, is_err=True)
        _msg("           " + " ".join(cmd), show_progress=show_progress, is_err=True)
        raise

def _is_nonempty_fasta(path: str) -> bool:
    try:
        if not os.path.isfile(path) or os.path.getsize(path) <= 0:
            return False
        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for line in fh:
                if line.startswith(">"):
                    return True
        return False
    except Exception:
        return False

def _classify_replicon(subdf: pd.DataFrame) -> Dict[str, Any]:
    hits = subdf["hit_gene_ref"].astype(str)
    relaxase_genes = sorted(set(hits[hits.str.startswith("T4SS_MOB")].tolist()))
    has_relaxase = len(relaxase_genes) > 0
    mpf_hits = hits[(hits.str.startswith("T4SS_")) & (~hits.str.startswith("T4SS_MOB"))]
    mpf_genes = sorted(set(mpf_hits.tolist()))
    tot_mpf = len(mpf_genes)
    has_t4cp = mpf_hits.str.contains("t4cp", case=False, na=False).any()
    has_virb4 = mpf_hits.str.contains("virb4", case=False, na=False).any()
    extra_mpf_genes = [g for g in mpf_genes if ("t4cp" not in g.lower()) and ("virb4" not in g.lower())]
    extra_mpf = len(extra_mpf_genes)
    mpf_types = set()
    for mfqn in subdf["model_fqn"].astype(str).unique():
        if "type" in mfqn:
            idx = mfqn.find("type")
            subtype = mfqn[idx + 4:].split("/")[0]
            mpf_types.add(subtype)
    mpf_type: Optional[str] = None
    if len(mpf_types) == 1:
        mpf_type = list(mpf_types)[0]
    elif len(mpf_types) > 1:
        mpf_type = ",".join(sorted(mpf_types))
    extra_threshold_map = {"FA": 2, "FATA": 2, "B": 3, "C": 3, "F": 3, "G": 3, "I": 3, "T": 3}
    extra_min = extra_threshold_map.get(mpf_type, 3)
    if not has_relaxase:
        mobility = "pNT"
    else:
        conj_condition = has_t4cp and has_virb4 and extra_mpf >= extra_min and tot_mpf >= 5
        mobility = "pCONJ" if conj_condition else "pMOB"
    return {"has_relaxase": has_relaxase, "has_t4cp": has_t4cp, "has_virb4": has_virb4, "tot_mpf": tot_mpf, "extra_mpf": extra_mpf, "mpf_type": mpf_type, "mobility": mobility, "relaxase_genes": ",".join(relaxase_genes), "mpf_genes": ",".join(mpf_genes), "extra_mpf_genes": ",".join(extra_mpf_genes)}

def _empty_row(genome_id: str) -> Dict[str, Any]:
    return {"genome": genome_id, "replicon": genome_id, "has_relaxase": False, "has_t4cp": False, "has_virb4": False, "tot_mpf": 0, "extra_mpf": 0, "mpf_type": None, "mobility": "pNT", "relaxase_genes": "", "mpf_genes": "", "extra_mpf_genes": ""}

def _parse_all_systems(tsv_path: str, genome_id: str) -> List[Dict[str, Any]]:
    if not os.path.exists(tsv_path):
        return [_empty_row(genome_id)]
    try:
        df = pd.read_csv(tsv_path, sep="\t", comment="#")
    except pd.errors.EmptyDataError:
        return [_empty_row(genome_id)]
    if df.empty:
        return [_empty_row(genome_id)]
    results: List[Dict[str, Any]] = []
    for rep, subdf in df.groupby("replicon"):
        info = _classify_replicon(subdf)
        row = {"genome": genome_id, "replicon": rep}
        row.update(info)
        results.append(row)
    return results

def standardize_mobility_output(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "genome" in out.columns and "PLASMID" not in out.columns:
        out = out.rename(columns={"genome": "PLASMID"})
    elif "PLASMID" not in out.columns:
        out.insert(0, "PLASMID", "")
    return out

def summarize_mobility_calls(input_fastas: List[str], mobility_df: pd.DataFrame) -> pd.DataFrame:
    plasmid_names = [os.path.basename(p) for p in input_fastas]
    if mobility_df is None or mobility_df.empty:
        return pd.DataFrame({"PLASMID": plasmid_names, "MOB": ["NA"] * len(plasmid_names)})
    records = []
    for plasmid in plasmid_names:
        stem = os.path.splitext(plasmid)[0]
        sub = mobility_df.loc[mobility_df["PLASMID"].isin([plasmid, stem])] if "PLASMID" in mobility_df.columns else pd.DataFrame()
        mob_call = "NA"
        if not sub.empty and "mobility" in sub.columns:
            if (sub["mobility"] == "pCONJ").any():
                mob_call = "pCONJ"
            elif (sub["mobility"] == "pMOB").any():
                mob_call = "pMOB"
            else:
                mob_call = "pNT"
        records.append({"PLASMID": plasmid, "MOB": mob_call})
    return pd.DataFrame(records)

def run_mobility(
    genome_fastas: Iterable[str],
    outdir: str,
    models_dir: str,
    prodigal_mode: str = "meta",
    prodigal_g: Optional[int] = None,
    prodigal_quiet: bool = False,
    macsy_evalue_search: Optional[float] = None,
    macsy_i_evalue_sel: Optional[float] = None,
    macsy_coverage_profile: Optional[float] = None,
    macsy_no_cut_ga: bool = False,
    macsy_force: bool = False,
    show_progress: bool = True,
) -> pd.DataFrame:
    genome_fastas = list(genome_fastas)
    os.makedirs(outdir, exist_ok=True)
    all_rows: List[Dict[str, Any]] = []
    for gpath in genome_fastas:
        gpath = os.path.abspath(gpath)
        if not os.path.isfile(gpath):
            _msg(f"[mobility] WARN: input FASTA not found, skip: {gpath}", show_progress=show_progress, is_err=True)
            continue
        g_dir = os.path.dirname(gpath)
        g_name = os.path.basename(gpath)
        base_no_ext = os.path.splitext(g_name)[0]
        protein_faa = os.path.join(g_dir, f"{base_no_ext}protein.fasta")
        prodigal_cmd = ["prodigal", "-a", protein_faa, "-p", prodigal_mode, "-i", gpath]
        if prodigal_g is not None:
            prodigal_cmd += ["-g", str(prodigal_g)]
        if prodigal_quiet:
            prodigal_cmd.append("-q")
        try:
            _run_cmd(prodigal_cmd, show_progress=show_progress)
        except Exception:
            if not _is_nonempty_fasta(protein_faa):
                all_rows.append(_empty_row(base_no_ext))
                continue
        if not _is_nonempty_fasta(protein_faa):
            all_rows.append(_empty_row(base_no_ext))
            continue
        msf_outdir = os.path.join(outdir, base_no_ext)
        os.makedirs(msf_outdir, exist_ok=True)
        macsy_cmd = ["macsyfinder", "--models", "CONJScan/Plasmids", "all", "--sequence-db", protein_faa, "--db-type", "ordered_replicon", "--replicon-topology", "circular", "--models-dir", models_dir, "-o", msf_outdir]
        if macsy_evalue_search is not None:
            macsy_cmd += ["--e-value-search", str(macsy_evalue_search)]
        if macsy_i_evalue_sel is not None:
            macsy_cmd += ["--i-evalue-sel", str(macsy_i_evalue_sel)]
        if macsy_coverage_profile is not None:
            macsy_cmd += ["--coverage-profile", str(macsy_coverage_profile)]
        if macsy_no_cut_ga:
            macsy_cmd.append("--no-cut-ga")
        if macsy_force:
            macsy_cmd.append("--force")
        _run_cmd(macsy_cmd, show_progress=show_progress)
        all_rows.extend(_parse_all_systems(os.path.join(msf_outdir, "all_systems.tsv"), base_no_ext))
    return standardize_mobility_output(pd.DataFrame(all_rows)) if all_rows else pd.DataFrame()
