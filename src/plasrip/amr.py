#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
AMR calling module: run abricate on plasmid FASTA files and merge results.
"""

import argparse
import os
import sys
import subprocess
from typing import List, Optional

import pandas as pd
from tqdm import tqdm

VERBOSE = True

def log(msg: str) -> None:
    if VERBOSE:
        sys.stderr.write(msg + "\n")

_ABRICATE_SUPPORTS_THREADS: Optional[bool] = None
_ABRICATE_THREADS_WARNED: bool = False

def _abricate_supports_threads() -> bool:
    global _ABRICATE_SUPPORTS_THREADS
    if _ABRICATE_SUPPORTS_THREADS is not None:
        return _ABRICATE_SUPPORTS_THREADS
    try:
        proc = subprocess.run(["abricate", "--help"], capture_output=True, text=True, check=False)
        text = (proc.stdout or "") + (proc.stderr or "")
        _ABRICATE_SUPPORTS_THREADS = "--threads" in text
    except Exception:
        _ABRICATE_SUPPORTS_THREADS = False
    return _ABRICATE_SUPPORTS_THREADS

def run_abricate_on_fasta(
    fasta: str,
    db: str,
    threads: int,
    minid: Optional[float] = None,
    mincov: Optional[float] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    global _ABRICATE_THREADS_WARNED
    cmd: List[str] = ["abricate", "--db", db]
    if threads and threads > 1 and _abricate_supports_threads():
        cmd += ["--threads", str(threads)]
    elif threads and threads > 1 and not _abricate_supports_threads():
        if not _ABRICATE_THREADS_WARNED and verbose:
            log("[WARN] abricate in this environment does not support '--threads'; running AMR step single-threaded.")
            _ABRICATE_THREADS_WARNED = True
    if minid is not None:
        cmd += ["--minid", str(minid)]
    if mincov is not None:
        cmd += ["--mincov", str(mincov)]
    cmd.append(fasta)
    if verbose:
        log(f"[RUN ] {' '.join(cmd)}")
    try:
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        if verbose:
            log(f"[ERROR] abricate failed on {os.path.basename(fasta)}")
            if e.stderr:
                log(e.stderr.strip())
        return pd.DataFrame()
    out = proc.stdout.strip()
    if not out:
        return pd.DataFrame()
    return pd.read_csv(pd.io.common.StringIO(out), sep="\t")

def summarize_amr_hits(input_fastas: List[str], amr_df: pd.DataFrame) -> pd.DataFrame:
    plasmid_names = [os.path.basename(p) for p in input_fastas]
    if amr_df is None or amr_df.empty:
        return pd.DataFrame({"PLASMID": plasmid_names, "AMR_hits": [0] * len(plasmid_names), "AMR_genes": [""] * len(plasmid_names)})
    gene_col = "GENE" if "GENE" in amr_df.columns else ("gene" if "gene" in amr_df.columns else None)
    records = []
    for plasmid in plasmid_names:
        sub = amr_df.loc[amr_df["PLASMID"] == plasmid] if "PLASMID" in amr_df.columns else pd.DataFrame()
        n_hits = int(sub.shape[0])
        genes = ";".join(pd.unique(sub[gene_col].astype(str))) if gene_col and not sub.empty else ""
        records.append({"PLASMID": plasmid, "AMR_hits": n_hits, "AMR_genes": genes})
    return pd.DataFrame(records)

def run_amr_pipeline(
    input_fastas: List[str],
    output: str,
    db: str = "ncbi",
    threads: int = 1,
    minid: Optional[float] = None,
    mincov: Optional[float] = None,
    show_progress: bool = True,
) -> pd.DataFrame:
    global VERBOSE
    VERBOSE = show_progress
    fastas: List[str] = []
    for p in input_fastas:
        if not os.path.isfile(p):
            log(f"[ERROR] FASTA not found: {p}")
            raise FileNotFoundError(p)
        fastas.append(os.path.abspath(p))
    if show_progress:
        log("========== abricate AMR calling ==========")
        log(f"[CFG ] db        : {db}")
        log(f"[CFG ] threads   : {threads}")
        if minid is not None:
            log(f"[CFG ] minid    : {minid}")
        if mincov is not None:
            log(f"[CFG ] mincov   : {mincov}")
        log(f"[CFG ] plasmids  : {len(fastas)} file(s)")
        log("==========================================")
    all_results: List[pd.DataFrame] = []
    iterator = tqdm(fastas, total=len(fastas), desc="Plasmids", unit="file", ncols=100) if show_progress else fastas
    for fa in iterator:
        base = os.path.basename(fa)
        df = run_abricate_on_fasta(fa, db, threads, minid=minid, mincov=mincov, verbose=show_progress)
        if df is None or df.empty:
            if show_progress:
                log(f"[HITS] {base}: 0 hit(s)")
            continue
        df = df.copy()
        df.insert(0, "PLASMID", base)
        all_results.append(df)
        if show_progress:
            log(f"[HITS] {base}: {len(df)} hit(s)")
    merged = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame(columns=["PLASMID"])
    merged.to_csv(output, sep="\t", index=False)
    if show_progress:
        log("==========================================")
        log(f"[OUT ] merged abricate hits: {output}")
        log(f"[OUT ] total hits          : {len(merged)}")
        log("==========================================")
        log("[DONE] abricate AMR calling finished")
        log("==========================================")
    return merged

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run abricate on plasmid FASTA files and merge results")
    parser.add_argument("-i", "--input", nargs="+", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-d", "--db", default="ncbi")
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("--minid", type=float, default=None)
    parser.add_argument("--mincov", type=float, default=None)
    return parser.parse_args()

def main() -> None:
    args = parse_args()
    run_amr_pipeline(args.input, args.output, db=args.db, threads=args.threads, minid=args.minid, mincov=args.mincov, show_progress=True)

if __name__ == "__main__":
    main()
