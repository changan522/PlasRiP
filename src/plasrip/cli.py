"""Top-level command line interface for PlasRiP v2 with unified terminal style."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from importlib.metadata import PackageNotFoundError, version

from .annotate import run_annotation
from .config import DEFAULT_SCORE_DELTA
from .core import run_core_analysis
from .pipeline import run_pipeline
from .predict import run_host_prediction
from .assign import run_ptu_assignment as _legacy_run_ptu_assignment
from .utils.console import HEADER_ART, print_header, print_info
from .utils.validate import check_dependencies, validate_model_dir, validate_rep_dir

try:
    __version__ = version("PlasRip")
except PackageNotFoundError:
    __version__ = "0.1.0"

def _json_default(obj):
    if isinstance(obj, Path):
        return str(obj)
    return str(obj)

def _cmd_predict(args) -> int:
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
    return 0

def _cmd_assign(args) -> int:
    _legacy_run_ptu_assignment(
        input_fastas=args.input,
        rep_dir=args.rep_dir,
        output=args.output,
        threshold=args.threshold,
        threads=args.threads,
        mash_cmd=args.mash,
        mash_k=args.mash_k,
        mash_s=args.mash_sketch_size,
        show_progress=not args.quiet,
    )
    return 0

def _cmd_core(args) -> int:
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
        mash_s=args.mash_sketch_size,
        output_all_scores=args.output_all_scores,
        cache_dir=args.cache_dir,
        keep_intermediate=args.keep_intermediate,
        threads=args.threads,
        overwrite=args.overwrite,
        log_file=args.log_file,
        quiet=args.quiet,
    )
    return 0

def _cmd_annotate(args) -> int:
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
    )
    return 0

def _cmd_run(args) -> int:
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
        mash_s=args.mash_sketch_size,
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
    return 0

def _cmd_version(args) -> int:
    _ = args
    print_header()
    print_info(f"PlasRiP {__version__}")
    print_info("Core modules: host-range prediction, PTU assignment")
    print_info("Annotation modules: AMR, replicon, mobility")
    return 0

def _cmd_models(args) -> int:
    report: dict[str, object] = {"plasrip_version": __version__, "dependencies": check_dependencies()}
    try:
        report["model_dir"] = validate_model_dir(args.model_dir)
    except Exception as exc:
        report["model_dir_error"] = str(exc)
    if args.rep_dir:
        try:
            report["rep_dir"] = validate_rep_dir(args.rep_dir)
        except Exception as exc:
            report["rep_dir_error"] = str(exc)
    if args.json:
        print(json.dumps(report, indent=2, ensure_ascii=False, default=_json_default))
    else:
        print_header()
        print_info(f"PlasRiP {__version__}")
        print_info("Dependencies:")
        for name, ok in report["dependencies"].items():
            print(f"  - {name}: {'OK' if ok else 'MISSING'}")
        if "model_dir" in report:
            print_info("Model directory:")
            for k, v in report["model_dir"].items():
                print(f"  - {k}: {v}")
        elif "model_dir_error" in report:
            print(f"Model directory error: {report['model_dir_error']}")
        if "rep_dir" in report:
            print_info(f"PTU sketch directory: {report['rep_dir']}")
        elif "rep_dir_error" in report:
            print(f"PTU sketch directory error: {report['rep_dir_error']}")
    return 0

def _add_common_predict_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input plasmid FASTA files")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Thread count")
    parser.add_argument("--topk", type=int, default=10, help="Top-k hosts to keep")
    parser.add_argument("--score-delta", type=float, default=DEFAULT_SCORE_DELTA, help="Keep all hosts with score >= top1_score - score_delta")
    parser.add_argument("--tax-level", choices=["species", "genus", "family"], default="species")
    parser.add_argument("--range-threshold-narrow", type=float, default=None)
    parser.add_argument("--range-threshold-broad", type=float, default=None)
    parser.add_argument("--output-all-scores", action="store_true")
    parser.add_argument("--cache-dir", default=None)
    parser.add_argument("--keep-intermediate", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--quiet", action="store_true")

def add_run_parser(subparsers) -> None:
    parser = subparsers.add_parser("run", help="Run full workflow")
    _add_common_predict_args(parser)
    parser.add_argument("--model-dir", required=False, default=None, help="Host-range model directory")
    parser.add_argument("--rep-dir", required=True, help="PTU member .msh directory")
    parser.add_argument("--ptu-threshold", type=float, default=0.05)
    parser.add_argument("--ptu-mash", default="mash")
    parser.add_argument("--mash-k", type=int, default=None)
    parser.add_argument("--mash-sketch-size", type=int, default=None)
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
    parser.set_defaults(func=_cmd_run)

def add_predict_parser(subparsers) -> None:
    parser = subparsers.add_parser("predict", help="Predict plasmid host range")
    _add_common_predict_args(parser)
    parser.add_argument("--model-dir", required=False, default=None, help="Host-range model directory")
    parser.set_defaults(func=_cmd_predict)

def add_assign_parser(subparsers) -> None:
    parser = subparsers.add_parser("assign", help="Assign plasmids to PTUs")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input plasmid FASTA files")
    parser.add_argument("-r", "--rep-dir", required=True, help="PTU member .msh directory")
    parser.add_argument("-o", "--output", required=True, help="PTU output TSV")
    parser.add_argument("--threshold", type=float, default=0.05, help="Mash distance cutoff")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Mash threads")
    parser.add_argument("--mash", default="mash", help="Mash executable")
    parser.add_argument("--mash-k", type=int, default=None, help="Mash sketch k-mer size")
    parser.add_argument("--mash-sketch-size", type=int, default=None, help="Mash sketch size")
    parser.add_argument("--quiet", action="store_true")
    parser.set_defaults(func=_cmd_assign)

def add_core_parser(subparsers) -> None:
    parser = subparsers.add_parser("core", help="Run host-range prediction + PTU assignment")
    _add_common_predict_args(parser)
    parser.add_argument("--model-dir", required=False, default=None, help="Host-range model directory")
    parser.add_argument("--rep-dir", required=True, help="PTU member .msh directory")
    parser.add_argument("--ptu-threshold", type=float, default=0.05)
    parser.add_argument("--ptu-mash", default="mash")
    parser.add_argument("--mash-k", type=int, default=None)
    parser.add_argument("--mash-sketch-size", type=int, default=None)
    parser.set_defaults(func=_cmd_core)

def add_annotate_parser(subparsers) -> None:
    parser = subparsers.add_parser("annotate", help="Run one-stop plasmid annotation")
    parser.add_argument("-i", "--input", nargs="+", required=True)
    parser.add_argument("-o", "--out-dir", required=True)
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
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--quiet", action="store_true")
    parser.set_defaults(func=_cmd_annotate)

def add_version_parser(subparsers) -> None:
    parser = subparsers.add_parser("version", help="Show version information")
    parser.set_defaults(func=_cmd_version)

def add_models_parser(subparsers) -> None:
    parser = subparsers.add_parser("models", help="Check model/data/dependency availability")
    parser.add_argument("--model-dir", default=None)
    parser.add_argument("--rep-dir", default=None)
    parser.add_argument("--json", action="store_true")
    parser.set_defaults(func=_cmd_models)

class PlasRipArgumentParser(argparse.ArgumentParser):
    def format_help(self) -> str:
        return HEADER_ART + super().format_help()

def build_parser() -> argparse.ArgumentParser:
    parser = PlasRipArgumentParser(prog="plasrip", description="PlasRiP v2 command line interface")
    subparsers = parser.add_subparsers(dest="command")
    add_run_parser(subparsers)
    add_predict_parser(subparsers)
    add_assign_parser(subparsers)
    add_core_parser(subparsers)
    add_annotate_parser(subparsers)
    add_version_parser(subparsers)
    add_models_parser(subparsers)
    return parser

def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    func = getattr(args, "func", None)
    if func is None:
        parser.print_help()
        return 0
    return int(func(args) or 0)
