"""Dual-core analysis entry point for PlasRiP v2.

Precision fix version:
- keep per-plasmid execution
- print host-range details before PTU assignment
- avoid repeated ASCII header when nested under `run`
- normalize PLASMID keys with Path(...).stem
- merge PTU_info.tsv into core_summary.tsv with old-style summary field names
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .assign import run_ptu_assignment as _legacy_run_ptu_assignment
from .config import (
    DEFAULT_CORE_LOG_NAME,
    DEFAULT_CORE_SUMMARY_OUT,
    DEFAULT_HOSTRANGE_OUT,
    DEFAULT_HOSTRANGE_SET_OUT,
    DEFAULT_HOSTRANGE_TAXONOMY_OUT,
    DEFAULT_HOSTRANGE_TOP_HITS_OUT,
    DEFAULT_PTU_HITS_SUFFIX,
    DEFAULT_PTU_OUT,
    get_ptu_info_path,
)
from .predict import run_host_prediction
from .schemas import (
    CORE_SUMMARY_COLUMNS,
    PTU_COLUMNS,
    PTU_HITS_COLUMNS,
    PTU_INFO_COLUMNS,
    empty_ptu_df,
    empty_ptu_hits_df,
    reorder_columns,
)
from .utils.console import (
    print_done,
    print_error,
    print_final_summary,
    print_header,
    print_info,
    print_plasmid_header,
    print_step,
    print_warn,
)
from .utils.io import append_tsv, ensure_dir, reset_tsv
from .utils.logging import close_logger, setup_logger


def _norm_plasmid_key(x: object) -> str:
    return Path(str(x)).stem if x is not None else ""


def _safe_read_ptu_info() -> pd.DataFrame:
    path = get_ptu_info_path()
    if not path.is_file():
        return pd.DataFrame(columns=PTU_INFO_COLUMNS)
    try:
        df = pd.read_csv(path, sep="\t")
    except Exception:
        return pd.DataFrame(columns=PTU_INFO_COLUMNS)
    if "PTU" not in df.columns:
        return pd.DataFrame(columns=PTU_INFO_COLUMNS)

    rename_map = {
        "n_plasmids": "PTU_n_plasmids",
        "host_range_level": "PTU_host_range_level",
        "host_species_range": "PTU_host_species_range",
        "mob_dominant": "PTU_mob_dominant",
        "avg_genome_size_bp": "PTU_avg_genome_size_bp",
        "genome_size_min_bp": "PTU_genome_size_min_bp",
        "genome_size_max_bp": "PTU_genome_size_max_bp",
        "genome_size_range_bp": "PTU_genome_size_range_bp",
        "avg_AMR_genes_per_plasmid": "PTU_avg_AMR_genes_per_plasmid",
        "unique_AMR_classes": "PTU_unique_AMR_classes",
    }
    df = df.rename(columns=rename_map)
    return reorder_columns(df.copy(), PTU_INFO_COLUMNS, keep_extra=True)


def _standardize_ptu_df(ptu_df: pd.DataFrame) -> pd.DataFrame:
    if ptu_df is None or ptu_df.empty:
        return empty_ptu_df()

    out = ptu_df.copy()
    rename_map = {
        "query": "PLASMID",
        "n_ptu_hits": "PTU_hits",
        "ptu_assigned": "PTU",
        "PTU_assign_method": "how_assigned",
        "PTU_closest_rep_dist": "closest_rep_dist",
    }
    out = out.rename(columns=rename_map)
    if "PLASMID" in out.columns:
        out["PLASMID"] = out["PLASMID"].map(_norm_plasmid_key)
    return reorder_columns(out, PTU_COLUMNS, keep_extra=True)


def _standardize_ptu_hits_df(ptu_hits_df: pd.DataFrame) -> pd.DataFrame:
    if ptu_hits_df is None or ptu_hits_df.empty:
        return empty_ptu_hits_df()
    out = ptu_hits_df.copy().rename(columns={"plasmid": "PLASMID"})
    if "PLASMID" in out.columns:
        out["PLASMID"] = out["PLASMID"].map(_norm_plasmid_key)
    return reorder_columns(out, PTU_HITS_COLUMNS, keep_extra=True)


def build_core_summary(hostrange_df: pd.DataFrame, ptu_df: pd.DataFrame) -> pd.DataFrame:
    left = hostrange_df.copy() if hostrange_df is not None else pd.DataFrame()
    right = ptu_df.copy() if ptu_df is not None else pd.DataFrame()

    if "PLASMID" in left.columns:
        left["PLASMID"] = left["PLASMID"].map(_norm_plasmid_key)
    if "PLASMID" in right.columns:
        right["PLASMID"] = right["PLASMID"].map(_norm_plasmid_key)

    if left.empty or "PLASMID" not in left.columns:
        if "PLASMID" in right.columns and not right.empty:
            left = pd.DataFrame({"PLASMID": right["PLASMID"].dropna().astype(str).drop_duplicates().tolist()})
        else:
            left = pd.DataFrame(columns=["PLASMID"])

    ptu_sub = right[[c for c in ["PLASMID", "PTU", "how_assigned", "PTU_hits", "closest_rep_dist"] if c in right.columns]].copy()
    if ptu_sub.empty or "PLASMID" not in ptu_sub.columns:
        summary = left.copy()
    else:
        summary = left.merge(ptu_sub, on="PLASMID", how="left")

    ptu_info = _safe_read_ptu_info()
    if not ptu_info.empty and "PTU" in summary.columns:
        summary = summary.merge(ptu_info, on="PTU", how="left")

    return reorder_columns(summary, CORE_SUMMARY_COLUMNS, keep_extra=False)


def _reset_core_outputs(out_root: Path, output_all_scores: bool = False) -> None:
    reset_tsv(out_root / DEFAULT_HOSTRANGE_OUT)
    reset_tsv(out_root / DEFAULT_HOSTRANGE_TOP_HITS_OUT)
    reset_tsv(out_root / DEFAULT_HOSTRANGE_TAXONOMY_OUT)
    reset_tsv(out_root / DEFAULT_HOSTRANGE_SET_OUT)
    reset_tsv(out_root / DEFAULT_PTU_OUT)
    reset_tsv(out_root / f"{DEFAULT_PTU_OUT}{DEFAULT_PTU_HITS_SUFFIX}")
    reset_tsv(out_root / DEFAULT_CORE_SUMMARY_OUT)
    if output_all_scores:
        from .config import DEFAULT_HOSTRANGE_ALL_SCORES_OUT
        reset_tsv(out_root / DEFAULT_HOSTRANGE_ALL_SCORES_OUT)


def _append_core_outputs(results: dict, out_root: Path, output_all_scores: bool = False) -> None:
    append_tsv(results.get("hostrange", pd.DataFrame()), out_root / DEFAULT_HOSTRANGE_OUT)
    append_tsv(results.get("top_hits", pd.DataFrame()), out_root / DEFAULT_HOSTRANGE_TOP_HITS_OUT)
    append_tsv(results.get("taxonomy", pd.DataFrame()), out_root / DEFAULT_HOSTRANGE_TAXONOMY_OUT)
    append_tsv(results.get("set_df", pd.DataFrame()), out_root / DEFAULT_HOSTRANGE_SET_OUT)
    append_tsv(results.get("ptu", pd.DataFrame()), out_root / DEFAULT_PTU_OUT)
    append_tsv(results.get("ptu_hits", pd.DataFrame()), out_root / f"{DEFAULT_PTU_OUT}{DEFAULT_PTU_HITS_SUFFIX}")
    append_tsv(results.get("core_summary", pd.DataFrame()), out_root / DEFAULT_CORE_SUMMARY_OUT)
    if output_all_scores and "all_scores" in results:
        from .config import DEFAULT_HOSTRANGE_ALL_SCORES_OUT
        append_tsv(results["all_scores"], out_root / DEFAULT_HOSTRANGE_ALL_SCORES_OUT)


def run_core_analysis(
    input_fastas: list[str],
    model_dir: str | None,
    rep_dir: str,
    out_dir: str,
    topk: int = 10,
    score_delta: float = 0.20,
    tax_level: str = "species",
    range_threshold_narrow: float | None = None,
    range_threshold_broad: float | None = None,
    ptu_threshold: float = 0.05,
    ptu_mash: str = "mash",
    mash_k: int | None = None,
    mash_s: int | None = None,
    output_all_scores: bool = False,
    cache_dir: str | None = None,
    keep_intermediate: bool = False,
    threads: int = 1,
    overwrite: bool = False,
    log_file: str | None = None,
    quiet: bool = False,
    show_header: bool = True,
) -> dict:
    out_root = ensure_dir(out_dir, overwrite=overwrite)
    logger = setup_logger("plasrip.core", log_file or str(out_root / DEFAULT_CORE_LOG_NAME), quiet=True)

    hostrange_frames: list[pd.DataFrame] = []
    top_hit_frames: list[pd.DataFrame] = []
    taxonomy_frames: list[pd.DataFrame] = []
    set_frames: list[pd.DataFrame] = []
    ptu_frames: list[pd.DataFrame] = []
    ptu_hit_frames: list[pd.DataFrame] = []
    core_frames: list[pd.DataFrame] = []
    all_score_frames: list[pd.DataFrame] = []
    total = len(input_fastas)
    completed = 0
    failed = 0

    try:
        _reset_core_outputs(out_root, output_all_scores=output_all_scores)
        if not quiet:
            if show_header:
                print_header()
            print_info("core now runs in per-plasmid mode")

        for idx, fasta in enumerate(input_fastas, start=1):
            name = Path(fasta).name
            workdir = out_root / "_per_plasmid" / Path(fasta).stem / "core"
            workdir.mkdir(parents=True, exist_ok=True)

            if not quiet:
                print_plasmid_header(idx, total, fasta)
            logger.info("Start plasmid: %s", fasta)

            try:
                if not quiet:
                    print_step(name, "Host-range prediction")
                host_results = run_host_prediction(
                    input_fastas=[fasta],
                    model_dir=model_dir,
                    out_dir=str(workdir / "predict"),
                    topk=topk,
                    score_delta=score_delta,
                    tax_level=tax_level,
                    range_threshold_narrow=range_threshold_narrow,
                    range_threshold_broad=range_threshold_broad,
                    output_all_scores=output_all_scores,
                    cache_dir=str(workdir / "cache" / "hostrange") if cache_dir else None,
                    keep_intermediate=keep_intermediate,
                    threads=threads,
                    overwrite=True,
                    log_file=None,
                    quiet=True,
                )
                hr = host_results.get("hostrange", pd.DataFrame())
                if not hr.empty and not quiet:
                    row = hr.iloc[0]
                    print_done(
                        f"{name} :: top1={row.get('Predicted_host_top1', 'NA')} "
                        f"set_size={row.get('Predicted_host_set_size', 'NA')} "
                        f"range={row.get('Predicted_host_range_class', 'NA')}"
                    )
                elif not quiet:
                    print_warn(f"{name} :: no host-range summary row produced")

                if not quiet:
                    print_step(name, "PTU assignment")
                ptu_out = workdir / "ptu.tsv"
                ptu_raw = _legacy_run_ptu_assignment(
                    input_fastas=[fasta],
                    rep_dir=rep_dir,
                    output=str(ptu_out),
                    threshold=ptu_threshold,
                    threads=threads,
                    mash_cmd=ptu_mash,
                    mash_k=mash_k,
                    mash_s=mash_s,
                    show_progress=False,
                )
                ptu_df = _standardize_ptu_df(ptu_raw)
                hits_path = Path(str(ptu_out) + DEFAULT_PTU_HITS_SUFFIX)
                if hits_path.is_file():
                    ptu_hits = pd.read_csv(hits_path, sep="\t")
                    ptu_hits = _standardize_ptu_hits_df(ptu_hits)
                else:
                    ptu_hits = empty_ptu_hits_df()

                if not quiet:
                    if not ptu_df.empty:
                        prow = ptu_df.iloc[0]
                        print_done(
                            f"{name} :: PTU={prow.get('PTU', 'NA')} "
                            f"method={prow.get('how_assigned', 'NA')} "
                            f"hits={prow.get('PTU_hits', 'NA')}"
                        )
                    else:
                        print_warn(f"{name} :: no PTU row produced")

                results = dict(host_results)
                results["ptu"] = ptu_df
                results["ptu_hits"] = ptu_hits
                results["core_summary"] = build_core_summary(results.get("hostrange", pd.DataFrame()), ptu_df)

                _append_core_outputs(results, out_root, output_all_scores=output_all_scores)

                hostrange_frames.append(results.get("hostrange", pd.DataFrame()))
                top_hit_frames.append(results.get("top_hits", pd.DataFrame()))
                taxonomy_frames.append(results.get("taxonomy", pd.DataFrame()))
                set_frames.append(results.get("set_df", pd.DataFrame()))
                ptu_frames.append(ptu_df)
                ptu_hit_frames.append(ptu_hits)
                core_frames.append(results["core_summary"])
                if output_all_scores and "all_scores" in results:
                    all_score_frames.append(results["all_scores"])

                completed += 1
                if not quiet:
                    print_done(f"{name} :: core summary appended")
            except Exception as exc:  # pragma: no cover
                failed += 1
                logger.exception("Failed plasmid: %s", fasta)
                if not quiet:
                    print_error(f"{name} :: {exc}")
                continue

        if not quiet:
            print_final_summary(total=total, completed=completed, failed=failed, out_dir=out_root)

        return {
            "hostrange": pd.concat(hostrange_frames, ignore_index=True) if hostrange_frames else pd.DataFrame(),
            "top_hits": pd.concat(top_hit_frames, ignore_index=True) if top_hit_frames else pd.DataFrame(),
            "taxonomy": pd.concat(taxonomy_frames, ignore_index=True) if taxonomy_frames else pd.DataFrame(),
            "set_df": pd.concat(set_frames, ignore_index=True) if set_frames else pd.DataFrame(),
            "ptu": pd.concat(ptu_frames, ignore_index=True) if ptu_frames else empty_ptu_df(),
            "ptu_hits": pd.concat(ptu_hit_frames, ignore_index=True) if ptu_hit_frames else empty_ptu_hits_df(),
            "core_summary": pd.concat(core_frames, ignore_index=True) if core_frames else pd.DataFrame(columns=CORE_SUMMARY_COLUMNS),
            "all_scores": pd.concat(all_score_frames, ignore_index=True) if all_score_frames else pd.DataFrame(),
        }
    finally:
        close_logger(logger)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run host-range prediction + PTU assignment")
    parser.add_argument("-i", "--input", nargs="+", required=True)
    parser.add_argument("--model-dir", required=False, default=None)
    parser.add_argument("--rep-dir", required=True)
    parser.add_argument("-o", "--out-dir", required=True)
    parser.add_argument("--topk", type=int, default=10)
    parser.add_argument("--score-delta", type=float, default=0.20)
    parser.add_argument("--tax-level", choices=["species", "genus", "family"], default="species")
    parser.add_argument("--range-threshold-narrow", type=float, default=None)
    parser.add_argument("--range-threshold-broad", type=float, default=None)
    parser.add_argument("--ptu-threshold", type=float, default=0.05)
    parser.add_argument("--ptu-mash", default="mash")
    parser.add_argument("--mash-k", type=int, default=None)
    parser.add_argument("--mash-sketch-size", dest="mash_s", type=int, default=None)
    parser.add_argument("--output-all-scores", action="store_true")
    parser.add_argument("--cache-dir", default=None)
    parser.add_argument("--keep-intermediate", action="store_true")
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--quiet", action="store_true")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    run_core_analysis(
        input_fastas=args.input,
        model_dir=args.model_dir,
        rep_dir=args.rep_dir,
        out_dir=args.out_dir,
        topk=args.topk,
        score_delta=args.score_delta,
        tax_level=args.tax_level,
        range_threshold_narrow=args.range_threshold_narrow,
        range_threshold_broad=args.range_threshold_broad,
        ptu_threshold=args.ptu_threshold,
        ptu_mash=args.ptu_mash,
        mash_k=args.mash_k,
        mash_s=args.mash_s,
        output_all_scores=args.output_all_scores,
        cache_dir=args.cache_dir,
        keep_intermediate=args.keep_intermediate,
        threads=args.threads,
        overwrite=args.overwrite,
        log_file=args.log_file,
        quiet=args.quiet,
        show_header=True,
    )


if __name__ == "__main__":
    main()
