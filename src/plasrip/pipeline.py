"""Full workflow orchestrator for PlasRiP v2.

Minimal separator fix:
- preserve the existing refined top-level content
- only add blank-line + horizontal-rule separators between major sections
- keep nested modules silent
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

from .annotate import run_annotation
from .config import (
    DEFAULT_AMR_OUT,
    DEFAULT_ANNOTATION_SUMMARY_OUT,
    DEFAULT_CORE_SUMMARY_OUT,
    DEFAULT_HOSTRANGE_OUT,
    DEFAULT_HOSTRANGE_SET_OUT,
    DEFAULT_HOSTRANGE_TAXONOMY_OUT,
    DEFAULT_HOSTRANGE_TOP_HITS_OUT,
    DEFAULT_LOG_NAME,
    DEFAULT_MOBILITY_OUT,
    DEFAULT_PTU_HITS_SUFFIX,
    DEFAULT_PTU_OUT,
    DEFAULT_REPLICON_OUT,
    DEFAULT_SUMMARY_OUT,
)
from .core import run_core_analysis
from .schemas import SUMMARY_COLUMNS, empty_summary_df, reorder_columns
from .utils.console import (
    print_done,
    print_error,
    print_final_summary,
    print_header,
    print_info,
    print_plasmid_header,
    print_step,
)
from .utils.io import append_tsv, ensure_dir, reset_tsv
from .utils.logging import close_logger, setup_logger


def _norm_plasmid_key(x: object) -> str:
    return Path(str(x)).stem if x is not None else ""


def _first_row(df: pd.DataFrame) -> dict:
    if df is None or df.empty:
        return {}
    return df.iloc[0].to_dict()


def _fmt(v: object, default: str = "NA") -> str:
    if v is None:
        return default
    try:
        if pd.isna(v):
            return default
    except Exception:
        pass
    return str(v)


def _shorten(text: object, max_len: int = 90) -> str:
    s = _fmt(text, "")
    if len(s) <= max_len:
        return s
    return s[: max_len - 3] + "..."


def _section_break() -> None:
    sys.stderr.write("\n" + "-" * 70 + "\n")
    sys.stderr.flush()


def build_final_summary(core_summary_df: pd.DataFrame, annotation_summary_df: pd.DataFrame) -> pd.DataFrame:
    left = core_summary_df.copy() if core_summary_df is not None else pd.DataFrame()
    right = annotation_summary_df.copy() if annotation_summary_df is not None else pd.DataFrame()

    if "PLASMID" in left.columns:
        left["PLASMID"] = left["PLASMID"].map(_norm_plasmid_key)
    if "PLASMID" in right.columns:
        right["PLASMID"] = right["PLASMID"].map(_norm_plasmid_key)

    if left.empty and right.empty:
        return empty_summary_df()

    summary = left.merge(right, on="PLASMID", how="left", suffixes=("", "_ann"))
    for col in list(summary.columns):
        if col.endswith("_ann"):
            base = col[:-4]
            if base in summary.columns:
                summary.drop(columns=[col], inplace=True)
            else:
                summary.rename(columns={col: base}, inplace=True)

    return reorder_columns(summary, SUMMARY_COLUMNS, keep_extra=False)


def _reset_pipeline_outputs(out_root: Path, output_all_scores: bool = False) -> None:
    for name in [
        DEFAULT_HOSTRANGE_OUT,
        DEFAULT_HOSTRANGE_TOP_HITS_OUT,
        DEFAULT_HOSTRANGE_TAXONOMY_OUT,
        DEFAULT_HOSTRANGE_SET_OUT,
        DEFAULT_PTU_OUT,
        f"{DEFAULT_PTU_OUT}{DEFAULT_PTU_HITS_SUFFIX}",
        DEFAULT_CORE_SUMMARY_OUT,
        DEFAULT_AMR_OUT,
        DEFAULT_REPLICON_OUT,
        DEFAULT_MOBILITY_OUT,
        DEFAULT_ANNOTATION_SUMMARY_OUT,
        DEFAULT_SUMMARY_OUT,
    ]:
        reset_tsv(out_root / name)
    if output_all_scores:
        from .config import DEFAULT_HOSTRANGE_ALL_SCORES_OUT
        reset_tsv(out_root / DEFAULT_HOSTRANGE_ALL_SCORES_OUT)


def _append_results(results: dict, out_root: Path, output_all_scores: bool = False) -> None:
    append_tsv(results.get("hostrange", pd.DataFrame()), out_root / DEFAULT_HOSTRANGE_OUT)
    append_tsv(results.get("top_hits", pd.DataFrame()), out_root / DEFAULT_HOSTRANGE_TOP_HITS_OUT)
    append_tsv(results.get("taxonomy", pd.DataFrame()), out_root / DEFAULT_HOSTRANGE_TAXONOMY_OUT)
    append_tsv(results.get("set_df", pd.DataFrame()), out_root / DEFAULT_HOSTRANGE_SET_OUT)
    append_tsv(results.get("ptu", pd.DataFrame()), out_root / DEFAULT_PTU_OUT)
    append_tsv(results.get("ptu_hits", pd.DataFrame()), out_root / f"{DEFAULT_PTU_OUT}{DEFAULT_PTU_HITS_SUFFIX}")
    append_tsv(results.get("core_summary", pd.DataFrame()), out_root / DEFAULT_CORE_SUMMARY_OUT)
    append_tsv(results.get("amr", pd.DataFrame()), out_root / DEFAULT_AMR_OUT)
    append_tsv(results.get("replicon", pd.DataFrame()), out_root / DEFAULT_REPLICON_OUT)
    append_tsv(results.get("mobility", pd.DataFrame()), out_root / DEFAULT_MOBILITY_OUT)
    append_tsv(results.get("annotation_summary", pd.DataFrame()), out_root / DEFAULT_ANNOTATION_SUMMARY_OUT)
    append_tsv(results.get("summary", pd.DataFrame()), out_root / DEFAULT_SUMMARY_OUT)
    if output_all_scores and "all_scores" in results:
        from .config import DEFAULT_HOSTRANGE_ALL_SCORES_OUT
        append_tsv(results["all_scores"], out_root / DEFAULT_HOSTRANGE_ALL_SCORES_OUT)


def run_pipeline(
    input_fastas: list[str],
    out_dir: str,
    model_dir: str | None,
    rep_dir: str,
    topk: int = 10,
    score_delta: float = 0.20,
    tax_level: str = "species",
    range_threshold_narrow: float | None = None,
    range_threshold_broad: float | None = None,
    ptu_threshold: float = 0.05,
    ptu_mash: str = "mash",
    mash_k: int | None = None,
    mash_s: int | None = None,
    with_amr: bool = True,
    with_replicon: bool = True,
    with_mobility: bool = True,
    amr_db: str = "ncbi",
    abricate_minid: float | None = None,
    abricate_mincov: float | None = None,
    prodigal_g: int | None = None,
    prodigal_quiet: bool = False,
    macsy_evalue_search: float | None = None,
    macsy_i_evalue_sel: float | None = None,
    macsy_coverage_profile: float | None = None,
    macsy_no_cut_ga: bool = False,
    macsy_force: bool = False,
    models_dir: str | None = None,
    output_all_scores: bool = False,
    cache_dir: str | None = None,
    keep_intermediate: bool = False,
    threads: int = 1,
    overwrite: bool = False,
    log_file: str | None = None,
    quiet: bool = False,
) -> dict:
    out_root = ensure_dir(out_dir, overwrite=overwrite)
    logger = setup_logger("plasrip.pipeline", log_file or str(out_root / DEFAULT_LOG_NAME), quiet=True)

    hostrange_frames: list[pd.DataFrame] = []
    top_hit_frames: list[pd.DataFrame] = []
    taxonomy_frames: list[pd.DataFrame] = []
    set_frames: list[pd.DataFrame] = []
    ptu_frames: list[pd.DataFrame] = []
    ptu_hit_frames: list[pd.DataFrame] = []
    core_frames: list[pd.DataFrame] = []
    amr_frames: list[pd.DataFrame] = []
    rep_frames: list[pd.DataFrame] = []
    mob_frames: list[pd.DataFrame] = []
    ann_frames: list[pd.DataFrame] = []
    summary_frames: list[pd.DataFrame] = []
    all_score_frames: list[pd.DataFrame] = []
    total = len(input_fastas)
    completed = 0
    failed = 0

    try:
        _reset_pipeline_outputs(out_root, output_all_scores=output_all_scores)
        if not quiet:
            print_header()
            print_info("PlasRiP run started")
            print_info(f"total plasmids : {total}")
            print_info("mode           : per-plasmid workflow")

        for idx, fasta in enumerate(input_fastas, start=1):
            name = Path(fasta).name
            workdir = out_root / "_per_plasmid" / Path(fasta).stem
            workdir.mkdir(parents=True, exist_ok=True)

            if not quiet:
                sys.stderr.write("\n")
                sys.stderr.flush()
                print_plasmid_header(idx, total, fasta)
            logger.info("Start plasmid: %s", fasta)

            try:
                core_results = run_core_analysis(
                    input_fastas=[fasta],
                    model_dir=model_dir,
                    rep_dir=rep_dir,
                    out_dir=str(workdir / "core"),
                    topk=topk,
                    score_delta=score_delta,
                    tax_level=tax_level,
                    range_threshold_narrow=range_threshold_narrow,
                    range_threshold_broad=range_threshold_broad,
                    ptu_threshold=ptu_threshold,
                    ptu_mash=ptu_mash,
                    mash_k=mash_k,
                    mash_s=mash_s,
                    output_all_scores=output_all_scores,
                    cache_dir=str(workdir / "cache" / "core") if cache_dir else None,
                    keep_intermediate=keep_intermediate,
                    threads=threads,
                    overwrite=True,
                    log_file=None,
                    quiet=True,
                    show_header=False,
                )
                core_row = _first_row(core_results.get("core_summary", pd.DataFrame()))

                if not quiet:
                    print_step(name, "Host-range prediction")
                    print_done(f"top1                  : {_fmt(core_row.get('Predicted_host_top1'))}")
                    print_done(f"host set size         : {_fmt(core_row.get('Predicted_host_set_size'))}")
                    print_done(f"range class           : {_fmt(core_row.get('Predicted_host_range_class'))}")
                    print_info(f"score delta           : {score_delta}")
                    _section_break()

                    print_step(name, "PTU assignment")
                    print_done(f"PTU                   : {_fmt(core_row.get('PTU'))}")
                    print_done(f"assignment method     : {_fmt(core_row.get('how_assigned'))}")
                    print_done(f"PTU hits              : {_fmt(core_row.get('PTU_hits'))}")
                    print_info(f"closest rep dist      : {_fmt(core_row.get('closest_rep_dist'))}")
                    _section_break()

                ann_results = run_annotation(
                    input_fastas=[fasta],
                    out_dir=str(workdir / "annotate"),
                    with_amr=with_amr,
                    with_replicon=with_replicon,
                    with_mobility=with_mobility,
                    amr_db=amr_db,
                    abricate_minid=abricate_minid,
                    abricate_mincov=abricate_mincov,
                    prodigal_g=prodigal_g,
                    prodigal_quiet=prodigal_quiet,
                    macsy_evalue_search=macsy_evalue_search,
                    macsy_i_evalue_sel=macsy_i_evalue_sel,
                    macsy_coverage_profile=macsy_coverage_profile,
                    macsy_no_cut_ga=macsy_no_cut_ga,
                    macsy_force=macsy_force,
                    models_dir=models_dir,
                    threads=threads,
                    overwrite=True,
                    log_file=None,
                    quiet=True,
                    show_header=False,
                )
                ann_row = _first_row(ann_results.get("annotation_summary", pd.DataFrame()))

                if not quiet:
                    print_step(name, "AMR detection")
                    print_done(f"AMR hits              : {_fmt(ann_row.get('AMR_hits'))}")
                    genes = _shorten(ann_row.get("AMR_genes"), max_len=90)
                    if genes:
                        print_info(f"AMR genes             : {genes}")
                    _section_break()

                    print_step(name, "Replicon typing")
                    print_done(f"replicon hits         : {_fmt(ann_row.get('REPLICON_hits'))}")
                    rep_genes = _shorten(ann_row.get("REPLICON_genes"), max_len=90)
                    if rep_genes:
                        print_info(f"replicon genes        : {rep_genes}")
                    _section_break()

                    print_step(name, "Mobility prediction")
                    print_done(f"mobility class        : {_fmt(ann_row.get('MOB'))}")
                    _section_break()

                results = {**core_results, **ann_results}
                results["summary"] = build_final_summary(
                    core_results.get("core_summary", pd.DataFrame()),
                    ann_results.get("annotation_summary", pd.DataFrame()),
                )

                _append_results(results, out_root, output_all_scores=output_all_scores)

                hostrange_frames.append(results.get("hostrange", pd.DataFrame()))
                top_hit_frames.append(results.get("top_hits", pd.DataFrame()))
                taxonomy_frames.append(results.get("taxonomy", pd.DataFrame()))
                set_frames.append(results.get("set_df", pd.DataFrame()))
                ptu_frames.append(results.get("ptu", pd.DataFrame()))
                ptu_hit_frames.append(results.get("ptu_hits", pd.DataFrame()))
                core_frames.append(results.get("core_summary", pd.DataFrame()))
                amr_frames.append(results.get("amr", pd.DataFrame()))
                rep_frames.append(results.get("replicon", pd.DataFrame()))
                mob_frames.append(results.get("mobility", pd.DataFrame()))
                ann_frames.append(results.get("annotation_summary", pd.DataFrame()))
                summary_frames.append(results.get("summary", pd.DataFrame()))
                if output_all_scores and "all_scores" in results:
                    all_score_frames.append(results["all_scores"])

                completed += 1
                if not quiet:
                    print_step(name, "Final summary")
                    print_done("core summary          : appended")
                    print_done("annotation summary    : appended")
                    print_done("final summary         : appended")
                    sys.stderr.write("\n")
                    sys.stderr.flush()
            except Exception as exc:
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
            "ptu": pd.concat(ptu_frames, ignore_index=True) if ptu_frames else pd.DataFrame(),
            "ptu_hits": pd.concat(ptu_hit_frames, ignore_index=True) if ptu_hit_frames else pd.DataFrame(),
            "core_summary": pd.concat(core_frames, ignore_index=True) if core_frames else pd.DataFrame(),
            "amr": pd.concat(amr_frames, ignore_index=True) if amr_frames else pd.DataFrame(),
            "replicon": pd.concat(rep_frames, ignore_index=True) if rep_frames else pd.DataFrame(),
            "mobility": pd.concat(mob_frames, ignore_index=True) if mob_frames else pd.DataFrame(),
            "annotation_summary": pd.concat(ann_frames, ignore_index=True) if ann_frames else pd.DataFrame(),
            "summary": pd.concat(summary_frames, ignore_index=True) if summary_frames else empty_summary_df(),
            "all_scores": pd.concat(all_score_frames, ignore_index=True) if all_score_frames else pd.DataFrame(),
        }
    finally:
        close_logger(logger)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the full PlasRiP v2 workflow")
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
    parser.add_argument("--with-amr", dest="with_amr", action="store_true", default=True)
    parser.add_argument("--no-amr", dest="with_amr", action="store_false")
    parser.add_argument("--with-replicon", dest="with_replicon", action="store_true", default=True)
    parser.add_argument("--no-replicon", dest="with_replicon", action="store_false")
    parser.add_argument("--with-mobility", dest="with_mobility", action="store_true", default=True)
    parser.add_argument("--no-mobility", dest="with_mobility", action="store_false")
    parser.add_argument("--amr-db", default="ncbi")
    parser.add_argument("--abricate-minid", type=float, default=None)
    parser.add_argument("--abricate-mincov", type=float, default=None)
    parser.add_argument("--prodigal-g", type=int, default=None)
    parser.add_argument("--prodigal-quiet", action="store_true")
    parser.add_argument("--macsy-evalue-search", type=float, default=None)
    parser.add_argument("--macsy-i-evalue-sel", type=float, default=None)
    parser.add_argument("--macsy-coverage-profile", type=float, default=None)
    parser.add_argument("--macsy-no-cut-ga", action="store_true")
    parser.add_argument("--macsy-force", action="store_true")
    parser.add_argument("--models-dir", default=None)
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
    run_pipeline(
        input_fastas=args.input,
        out_dir=args.out_dir,
        model_dir=args.model_dir,
        rep_dir=args.rep_dir,
        topk=args.topk,
        score_delta=args.score_delta,
        tax_level=args.tax_level,
        range_threshold_narrow=args.range_threshold_narrow,
        range_threshold_broad=args.range_threshold_broad,
        ptu_threshold=args.ptu_threshold,
        ptu_mash=args.ptu_mash,
        mash_k=args.mash_k,
        mash_s=args.mash_s,
        with_amr=args.with_amr,
        with_replicon=args.with_replicon,
        with_mobility=args.with_mobility,
        amr_db=args.amr_db,
        abricate_minid=args.abricate_minid,
        abricate_mincov=args.abricate_mincov,
        prodigal_g=args.prodigal_g,
        prodigal_quiet=args.prodigal_quiet,
        macsy_evalue_search=args.macsy_evalue_search,
        macsy_i_evalue_sel=args.macsy_i_evalue_sel,
        macsy_coverage_profile=args.macsy_coverage_profile,
        macsy_no_cut_ga=args.macsy_no_cut_ga,
        macsy_force=args.macsy_force,
        models_dir=args.models_dir,
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
