"""High-level host-range prediction entry point for PlasRiP v2 (3.2 unified)."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .config import (
    DEFAULT_HOSTRANGE_ALL_SCORES_OUT,
    DEFAULT_HOSTRANGE_OUT,
    DEFAULT_HOSTRANGE_SET_OUT,
    DEFAULT_HOSTRANGE_TAXONOMY_OUT,
    DEFAULT_HOSTRANGE_TOP_HITS_OUT,
    DEFAULT_MODEL_FILENAME,
    DEFAULT_PREDICT_LOG_NAME,
    DEFAULT_SCORE_DELTA,
    DEFAULT_TAX_LEVEL,
    DEFAULT_TOPK,
)
from .hostrange_features import extract_plasmid_features
from .hostrange_model import score_hosts
from .hostrange_postprocess import postprocess_host_scores
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
from .utils.validate import validate_model_dir


def _prepare_predict_environment(model_dir: str | None, out_dir: str, cache_dir: str | None, overwrite: bool):
    out_root = ensure_dir(out_dir, overwrite=overwrite)
    model_files = validate_model_dir(model_dir)
    cache_root = Path(cache_dir).resolve() if cache_dir else out_root / "cache"
    cache_root.mkdir(parents=True, exist_ok=True)
    return out_root, cache_root, model_files


def _reset_predict_outputs(out_dir: Path, output_all_scores: bool = False) -> None:
    reset_tsv(out_dir / DEFAULT_HOSTRANGE_OUT)
    reset_tsv(out_dir / DEFAULT_HOSTRANGE_TOP_HITS_OUT)
    reset_tsv(out_dir / DEFAULT_HOSTRANGE_TAXONOMY_OUT)
    reset_tsv(out_dir / DEFAULT_HOSTRANGE_SET_OUT)
    if output_all_scores:
        reset_tsv(out_dir / DEFAULT_HOSTRANGE_ALL_SCORES_OUT)


def _append_predict_outputs(results: dict, out_dir: Path, output_all_scores: bool = False) -> None:
    append_tsv(results.get("hostrange", pd.DataFrame()), out_dir / DEFAULT_HOSTRANGE_OUT)
    append_tsv(results.get("top_hits", pd.DataFrame()), out_dir / DEFAULT_HOSTRANGE_TOP_HITS_OUT)
    append_tsv(results.get("taxonomy", pd.DataFrame()), out_dir / DEFAULT_HOSTRANGE_TAXONOMY_OUT)
    append_tsv(results.get("set_df", pd.DataFrame()), out_dir / DEFAULT_HOSTRANGE_SET_OUT)
    if output_all_scores and "all_scores" in results:
        append_tsv(results["all_scores"], out_dir / DEFAULT_HOSTRANGE_ALL_SCORES_OUT)


def run_host_prediction(
    input_fastas: list[str],
    model_dir: str | None,
    out_dir: str,
    topk: int = DEFAULT_TOPK,
    score_delta: float = DEFAULT_SCORE_DELTA,
    tax_level: str = DEFAULT_TAX_LEVEL,
    range_threshold_narrow: float | None = None,
    range_threshold_broad: float | None = None,
    output_all_scores: bool = False,
    cache_dir: str | None = None,
    keep_intermediate: bool = False,
    threads: int = 1,
    overwrite: bool = False,
    log_file: str | None = None,
    quiet: bool = False,
) -> dict:
    out_root, cache_root, model_files = _prepare_predict_environment(model_dir, out_dir, cache_dir, overwrite)
    model_root = str(model_files["root"])
    model_name = str(model_files.get("model_name") or DEFAULT_MODEL_FILENAME)
    logger = setup_logger("plasrip.predict", log_file or str(out_root / DEFAULT_PREDICT_LOG_NAME), quiet=True)

    hostrange_frames: list[pd.DataFrame] = []
    top_hit_frames: list[pd.DataFrame] = []
    taxonomy_frames: list[pd.DataFrame] = []
    set_frames: list[pd.DataFrame] = []
    all_score_frames: list[pd.DataFrame] = []
    total = len(input_fastas)
    completed = 0
    failed = 0

    try:
        _reset_predict_outputs(out_root, output_all_scores=output_all_scores)
        if not quiet:
            print_header()
            print_info("predict now runs in per-plasmid mode")
            print_info(f"host-range model (3.2): {model_name}")
            print_info(f"model directory       : {model_root}")

        for idx, fasta in enumerate(input_fastas, start=1):
            plasmid_key = Path(fasta).stem
            plasmid_cache = cache_root / plasmid_key
            plasmid_cache.mkdir(parents=True, exist_ok=True)
            if not quiet:
                print_plasmid_header(idx, total, fasta)
            logger.info("Start plasmid: %s", fasta)
            logger.info("Using host-range model: %s", model_name)
            try:
                if not quiet:
                    print_step(Path(fasta).name, "Host-range prediction")
                features = extract_plasmid_features(
                    input_fastas=[fasta],
                    model_dir=model_root,
                    cache_dir=str(plasmid_cache),
                    threads=threads,
                    keep_intermediate=keep_intermediate,
                )
                score_df = score_hosts(features=features, model_dir=model_root, threads=threads)
                results = postprocess_host_scores(
                    score_df=score_df,
                    model_dir=model_root,
                    topk=topk,
                    score_delta=score_delta,
                    tax_level=tax_level,
                    range_threshold_narrow=range_threshold_narrow,
                    range_threshold_broad=range_threshold_broad,
                )
                if output_all_scores:
                    results["all_scores"] = score_df

                _append_predict_outputs(results, out_root, output_all_scores=output_all_scores)

                if not results["hostrange"].empty:
                    top1 = str(results["hostrange"].iloc[0].get("Predicted_host_top1", ""))
                    if not quiet:
                        print_done(f"{Path(fasta).name} :: top1 host = {top1}")
                else:
                    if not quiet:
                        print_warn(f"{Path(fasta).name} :: no host prediction rows produced")

                hostrange_frames.append(results.get("hostrange", pd.DataFrame()))
                top_hit_frames.append(results.get("top_hits", pd.DataFrame()))
                taxonomy_frames.append(results.get("taxonomy", pd.DataFrame()))
                set_frames.append(results.get("set_df", pd.DataFrame()))
                if output_all_scores and "all_scores" in results:
                    all_score_frames.append(results["all_scores"])
                completed += 1
                logger.info("Completed plasmid: %s", fasta)
            except Exception as exc:  # pragma: no cover
                failed += 1
                logger.exception("Failed plasmid: %s", fasta)
                if not quiet:
                    print_error(f"{Path(fasta).name} :: {exc}")
                continue

        if not quiet:
            print_final_summary(total=total, completed=completed, failed=failed, out_dir=out_root)

        return {
            "hostrange": pd.concat(hostrange_frames, ignore_index=True) if hostrange_frames else pd.DataFrame(),
            "top_hits": pd.concat(top_hit_frames, ignore_index=True) if top_hit_frames else pd.DataFrame(),
            "taxonomy": pd.concat(taxonomy_frames, ignore_index=True) if taxonomy_frames else pd.DataFrame(),
            "set_df": pd.concat(set_frames, ignore_index=True) if set_frames else pd.DataFrame(),
            "all_scores": pd.concat(all_score_frames, ignore_index=True) if all_score_frames else pd.DataFrame(),
        }
    finally:
        close_logger(logger)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Predict plasmid host range with PlasRiP v2 (3.2 unified)")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input plasmid FASTA files")
    parser.add_argument("--model-dir", required=False, default=None, help="3.2 host-range model directory")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
    parser.add_argument("--topk", type=int, default=DEFAULT_TOPK)
    parser.add_argument("--score-delta", type=float, default=DEFAULT_SCORE_DELTA, help="Keep all hosts with score >= top1_score - score_delta")
    parser.add_argument("--tax-level", default=DEFAULT_TAX_LEVEL, choices=["species", "genus", "family"])
    parser.add_argument("--range-threshold-narrow", type=float, default=None)
    parser.add_argument("--range-threshold-broad", type=float, default=None)
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
    run_host_prediction(
        input_fastas=args.input,
        model_dir=args.model_dir,
        out_dir=args.out_dir,
        topk=args.topk,
        score_delta=args.score_delta,
        tax_level=args.tax_level,
        range_threshold_narrow=args.range_threshold_narrow,
        range_threshold_broad=args.range_threshold_broad,
        output_all_scores=args.output_all_scores,
        cache_dir=args.cache_dir,
        keep_intermediate=args.keep_intermediate,
        threads=args.threads,
        overwrite=args.overwrite,
        log_file=args.log_file,
        quiet=args.quiet,
    )


if __name__ == "__main__":
    main()
