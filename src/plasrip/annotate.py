"""High-level annotation entry point for PlasRiP v2.

Precision fix version:
- keep per-plasmid execution
- avoid repeated ASCII header when nested under `run`
- normalize annotation summary PLASMID to Path(...).stem for stable summary merging
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from .amr import run_amr_pipeline, summarize_amr_hits
from .config import (
    DEFAULT_AMR_DB,
    DEFAULT_AMR_OUT,
    DEFAULT_ANNOTATE_LOG_NAME,
    DEFAULT_ANNOTATION_SUMMARY_OUT,
    DEFAULT_MOBILITY_OUT,
    DEFAULT_REPLICON_OUT,
    get_conjscan_dir,
)
from .mobility import run_mobility, summarize_mobility_calls
from .replicon import run_replicon_typing, summarize_replicon_hits
from .schemas import (
    ANNOTATION_SUMMARY_COLUMNS,
    empty_annotation_summary_df,
    ensure_columns,
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
)
from .utils.fasta import validate_fastas
from .utils.io import append_tsv, ensure_dir, reset_tsv
from .utils.logging import close_logger, setup_logger
from .utils.validate import check_dependencies


def _norm_plasmid_key(x: object) -> str:
    return Path(str(x)).stem if x is not None else ""


def build_annotation_summary(
    input_fastas: list[str],
    amr_df: pd.DataFrame | None = None,
    replicon_df: pd.DataFrame | None = None,
    mobility_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    plasmids = [_norm_plasmid_key(p) for p in input_fastas]
    summary = pd.DataFrame({"PLASMID": plasmids})

    if amr_df is not None:
        amr_sum = summarize_amr_hits(input_fastas, amr_df)
        if "PLASMID" in amr_sum.columns:
            amr_sum["PLASMID"] = amr_sum["PLASMID"].map(_norm_plasmid_key)
        summary = summary.merge(amr_sum, on="PLASMID", how="left")
    else:
        summary["AMR_hits"] = 0
        summary["AMR_genes"] = ""

    if replicon_df is not None:
        rep_sum = summarize_replicon_hits(input_fastas, replicon_df)
        if "PLASMID" in rep_sum.columns:
            rep_sum["PLASMID"] = rep_sum["PLASMID"].map(_norm_plasmid_key)
        summary = summary.merge(rep_sum, on="PLASMID", how="left")
    else:
        summary["REPLICON_hits"] = 0
        summary["REPLICON_genes"] = ""

    if mobility_df is not None:
        mob_sum = summarize_mobility_calls(input_fastas, mobility_df)
        if "PLASMID" in mob_sum.columns:
            mob_sum["PLASMID"] = mob_sum["PLASMID"].map(_norm_plasmid_key)
        summary = summary.merge(mob_sum, on="PLASMID", how="left")
    else:
        summary["MOB"] = "NA"

    summary = ensure_columns(summary, ANNOTATION_SUMMARY_COLUMNS, fill_value="")
    return reorder_columns(summary, ANNOTATION_SUMMARY_COLUMNS, keep_extra=False)


def _reset_annotation_outputs(out_root: Path) -> None:
    reset_tsv(out_root / DEFAULT_AMR_OUT)
    reset_tsv(out_root / DEFAULT_REPLICON_OUT)
    reset_tsv(out_root / DEFAULT_MOBILITY_OUT)
    reset_tsv(out_root / DEFAULT_ANNOTATION_SUMMARY_OUT)


def _append_annotation_outputs(results: dict[str, pd.DataFrame], out_root: Path) -> None:
    if "amr" in results:
        append_tsv(results["amr"], out_root / DEFAULT_AMR_OUT)
    if "replicon" in results:
        append_tsv(results["replicon"], out_root / DEFAULT_REPLICON_OUT)
    if "mobility" in results:
        append_tsv(results["mobility"], out_root / DEFAULT_MOBILITY_OUT)
    if "annotation_summary" in results:
        append_tsv(results["annotation_summary"], out_root / DEFAULT_ANNOTATION_SUMMARY_OUT)


def run_annotation(
    input_fastas: list[str],
    out_dir: str,
    with_amr: bool = True,
    with_replicon: bool = True,
    with_mobility: bool = True,
    amr_db: str = DEFAULT_AMR_DB,
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
    threads: int = 1,
    overwrite: bool = False,
    log_file: str | None = None,
    quiet: bool = False,
    show_header: bool = True,
) -> dict[str, pd.DataFrame]:
    fastas = validate_fastas(input_fastas)
    out_root = ensure_dir(out_dir, overwrite=overwrite)
    logger = setup_logger("plasrip.annotate", log_file or str(out_root / DEFAULT_ANNOTATE_LOG_NAME), quiet=True)

    amr_frames: list[pd.DataFrame] = []
    rep_frames: list[pd.DataFrame] = []
    mob_frames: list[pd.DataFrame] = []
    ann_frames: list[pd.DataFrame] = []
    total = len(fastas)
    completed = 0
    failed = 0

    try:
        check_dependencies(
            need_abricate=(with_amr or with_replicon),
            need_prodigal=with_mobility,
            need_macsyfinder=with_mobility,
        )
        _reset_annotation_outputs(out_root)
        if not quiet:
            if show_header:
                print_header()
            print_info("annotate now runs in per-plasmid mode")

        for idx, fasta in enumerate(fastas, start=1):
            name = Path(fasta).name
            workdir = out_root / "_per_plasmid" / Path(fasta).stem / "annotate"
            workdir.mkdir(parents=True, exist_ok=True)

            if not quiet:
                print_plasmid_header(idx, total, fasta)
            logger.info("Start plasmid: %s", fasta)

            try:
                results: dict[str, pd.DataFrame] = {}

                if with_amr:
                    if not quiet:
                        print_step(name, "AMR detection")
                    results["amr"] = run_amr_pipeline(
                        input_fastas=[fasta],
                        output=str(workdir / DEFAULT_AMR_OUT),
                        db=amr_db,
                        threads=threads,
                        minid=abricate_minid,
                        mincov=abricate_mincov,
                        show_progress=False,
                    )
                    if not quiet:
                        print_done(f"{name} :: AMR done")

                if with_replicon:
                    if not quiet:
                        print_step(name, "Replicon typing")
                    results["replicon"] = run_replicon_typing(
                        input_fastas=[fasta],
                        output=str(workdir / DEFAULT_REPLICON_OUT),
                        threads=threads,
                        minid=abricate_minid,
                        mincov=abricate_mincov,
                        overwrite=True,
                        show_progress=False,
                    )
                    if not quiet:
                        print_done(f"{name} :: replicon done")

                if with_mobility:
                    if not quiet:
                        print_step(name, "Mobility prediction")
                    results["mobility"] = run_mobility(
                        genome_fastas=[fasta],
                        outdir=str(workdir / "mobility"),
                        models_dir=models_dir or str(get_conjscan_dir()),
                        prodigal_mode="meta",
                        prodigal_g=prodigal_g,
                        prodigal_quiet=prodigal_quiet,
                        macsy_evalue_search=macsy_evalue_search,
                        macsy_i_evalue_sel=macsy_i_evalue_sel,
                        macsy_coverage_profile=macsy_coverage_profile,
                        macsy_no_cut_ga=macsy_no_cut_ga,
                        macsy_force=macsy_force,
                        show_progress=False,
                    )
                    if not quiet:
                        mob_df = results["mobility"]
                        mob_call = "NA"
                        if mob_df is not None and not mob_df.empty and "mobility" in mob_df.columns:
                            if (mob_df["mobility"] == "pCONJ").any():
                                mob_call = "pCONJ"
                            elif (mob_df["mobility"] == "pMOB").any():
                                mob_call = "pMOB"
                            else:
                                mob_call = "pNT"
                        print_done(f"{name} :: mobility={mob_call}")

                if not results:
                    results["annotation_summary"] = empty_annotation_summary_df()
                else:
                    results["annotation_summary"] = build_annotation_summary(
                        input_fastas=[fasta],
                        amr_df=results.get("amr"),
                        replicon_df=results.get("replicon"),
                        mobility_df=results.get("mobility"),
                    )

                _append_annotation_outputs(results, out_root)

                if "amr" in results:
                    amr_frames.append(results["amr"])
                if "replicon" in results:
                    rep_frames.append(results["replicon"])
                if "mobility" in results:
                    mob_frames.append(results["mobility"])
                ann_frames.append(results["annotation_summary"])

                completed += 1
                if not quiet:
                    print_done(f"{name} :: annotation summary appended")
            except Exception as exc:  # pragma: no cover
                failed += 1
                logger.exception("Failed plasmid: %s", fasta)
                if not quiet:
                    print_error(f"{name} :: {exc}")
                continue

        if not quiet:
            print_final_summary(total=total, completed=completed, failed=failed, out_dir=out_root)

        return {
            "amr": pd.concat(amr_frames, ignore_index=True) if amr_frames else pd.DataFrame(),
            "replicon": pd.concat(rep_frames, ignore_index=True) if rep_frames else pd.DataFrame(),
            "mobility": pd.concat(mob_frames, ignore_index=True) if mob_frames else pd.DataFrame(),
            "annotation_summary": pd.concat(ann_frames, ignore_index=True) if ann_frames else empty_annotation_summary_df(),
        }
    finally:
        close_logger(logger)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run one-stop plasmid annotation for PlasRiP v2")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input plasmid FASTA files")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
    parser.add_argument("--with-amr", dest="with_amr", action="store_true", default=True, help="Run AMR detection")
    parser.add_argument("--no-amr", dest="with_amr", action="store_false", help="Disable AMR detection")
    parser.add_argument("--with-replicon", dest="with_replicon", action="store_true", default=True, help="Run replicon typing")
    parser.add_argument("--no-replicon", dest="with_replicon", action="store_false", help="Disable replicon typing")
    parser.add_argument("--with-mobility", dest="with_mobility", action="store_true", default=True, help="Run mobility prediction")
    parser.add_argument("--no-mobility", dest="with_mobility", action="store_false", help="Disable mobility prediction")
    parser.add_argument("--amr-db", default=DEFAULT_AMR_DB, help="Abricate AMR database name")
    parser.add_argument("--abricate-minid", type=float, default=None, help="Abricate minimum identity")
    parser.add_argument("--abricate-mincov", type=float, default=None, help="Abricate minimum coverage")
    parser.add_argument("--prodigal-g", type=int, default=None, help="Prodigal translation table")
    parser.add_argument("--prodigal-quiet", action="store_true", help="Run prodigal quietly")
    parser.add_argument("--macsy-evalue-search", type=float, default=None)
    parser.add_argument("--macsy-i-evalue-sel", type=float, default=None)
    parser.add_argument("--macsy-coverage-profile", type=float, default=None)
    parser.add_argument("--macsy-no-cut-ga", action="store_true")
    parser.add_argument("--macsy-force", action="store_true")
    parser.add_argument("--models-dir", default=None, help="MacSyFinder models directory")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Thread count")
    parser.add_argument("--overwrite", action="store_true", help="Reuse a non-empty output directory")
    parser.add_argument("--log-file", default=None, help="Optional log file path")
    parser.add_argument("--quiet", action="store_true", help="Suppress console progress output")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    run_annotation(
        input_fastas=args.input,
        out_dir=args.out_dir,
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
        threads=args.threads,
        overwrite=args.overwrite,
        log_file=args.log_file,
        quiet=args.quiet,
        show_header=True,
    )


if __name__ == "__main__":
    main()
