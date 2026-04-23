#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Replicon typing module for PlasRiP v2."""

from __future__ import annotations

import argparse
import os
from typing import List, Optional

import pandas as pd

from .amr import run_abricate_on_fasta

def standardize_replicon_output(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    if "PLASMID" not in out.columns:
        out.insert(0, "PLASMID", "")
    return out

def summarize_replicon_hits(input_fastas: List[str], replicon_df: pd.DataFrame) -> pd.DataFrame:
    plasmid_names = [os.path.basename(p) for p in input_fastas]
    if replicon_df is None or replicon_df.empty:
        return pd.DataFrame({"PLASMID": plasmid_names, "REPLICON_hits": [0] * len(plasmid_names), "REPLICON_genes": [""] * len(plasmid_names)})
    gene_col = "GENE" if "GENE" in replicon_df.columns else ("gene" if "gene" in replicon_df.columns else None)
    records = []
    for plasmid in plasmid_names:
        sub = replicon_df.loc[replicon_df["PLASMID"] == plasmid]
        n_hits = int(sub.shape[0])
        genes = ";".join(pd.unique(sub[gene_col].astype(str))) if gene_col and not sub.empty else ""
        records.append({"PLASMID": plasmid, "REPLICON_hits": n_hits, "REPLICON_genes": genes})
    return pd.DataFrame(records)

def run_plasmidfinder_on_fasta(
    fasta: str,
    threads: int = 1,
    minid: Optional[float] = None,
    mincov: Optional[float] = None,
    verbose: bool = True,
) -> pd.DataFrame:
    df = run_abricate_on_fasta(fasta=fasta, db="plasmidfinder", threads=threads, minid=minid, mincov=mincov, verbose=verbose)
    if df is None or df.empty:
        return pd.DataFrame()
    df = df.copy()
    df.insert(0, "PLASMID", os.path.basename(fasta))
    return standardize_replicon_output(df)

def run_replicon_typing(
    input_fastas: List[str],
    output: str,
    threads: int = 1,
    minid: Optional[float] = None,
    mincov: Optional[float] = None,
    overwrite: bool = False,
    show_progress: bool = True,
) -> pd.DataFrame:
    if os.path.exists(output) and not overwrite:
        raise FileExistsError(f"Output already exists: {output}. Use overwrite=True to replace it.")
    valid_fastas = [os.path.abspath(p) for p in input_fastas if os.path.isfile(p)]
    if len(valid_fastas) != len(input_fastas):
        missing = [p for p in input_fastas if not os.path.isfile(p)]
        raise FileNotFoundError(missing[0])
    all_results: List[pd.DataFrame] = []
    for fasta in valid_fastas:
        df = run_plasmidfinder_on_fasta(fasta=fasta, threads=threads, minid=minid, mincov=mincov, verbose=show_progress)
        if not df.empty:
            all_results.append(df)
    merged = pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame(columns=["PLASMID"])
    merged.to_csv(output, sep="\t", index=False)
    return merged

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run plasmidfinder-based replicon typing and append results per plasmid")
    parser.add_argument("-i", "--input", nargs="+", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("--minid", type=float, default=None)
    parser.add_argument("--mincov", type=float, default=None)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()

def main() -> None:
    args = parse_args()
    run_replicon_typing(args.input, args.output, threads=args.threads, minid=args.minid, mincov=args.mincov, overwrite=args.overwrite, show_progress=True)

if __name__ == "__main__":
    main()
