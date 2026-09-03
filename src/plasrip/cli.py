"""Top-level command line interface for PlasRiP."""

from __future__ import annotations

import argparse
import json
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from .annotate import run_annotation
from .assign import run_ptu_assignment
from .config import DEFAULT_SCORE_DELTA
from .core import run_core_analysis
from .database import validate_database_dir
from .exceptions import PlasRipError
from .pipeline import run_pipeline
from .predict import run_host_prediction
from .utils.console import HEADER_ART, print_error, print_header, print_info
from .utils.validate import check_dependencies, validate_model_dir

try:
    __version__ = version("PlasRip")
except PackageNotFoundError:
    __version__ = "0.2.1"


def _json_default(obj):
    if isinstance(obj, Path):
        return str(obj)
    return str(obj)


def _cmd_predict(args) -> int:
    run_host_prediction(
        input_fastas=args.input, model_dir=args.model_dir, out_dir=args.out_dir,
        topk=args.topk, score_delta=args.score_delta, tax_level=args.tax_level,
        range_threshold_narrow=args.range_threshold_narrow,
        range_threshold_broad=args.range_threshold_broad,
        output_all_scores=args.output_all_scores, cache_dir=args.cache_dir,
        keep_intermediate=args.keep_intermediate, threads=args.threads,
        overwrite=args.overwrite, log_file=args.log_file, quiet=args.quiet,
    )
    return 0


def _cmd_assign(args) -> int:
    run_ptu_assignment(
        input_fastas=args.input, database_dir=args.database, output=args.output,
        threshold=args.ptu_threshold, threads=args.threads, mash_cmd=args.mash,
        mash_k=args.mash_k, mash_s=args.mash_sketch_size,
        show_progress=not args.quiet,
    )
    return 0


def _cmd_core(args) -> int:
    run_core_analysis(
        input_fastas=args.input, model_dir=args.model_dir, database_dir=args.database,
        out_dir=args.out_dir, topk=args.topk, score_delta=args.score_delta,
        tax_level=args.tax_level, range_threshold_narrow=args.range_threshold_narrow,
        range_threshold_broad=args.range_threshold_broad,
        ptu_threshold=args.ptu_threshold, ptu_mash=args.ptu_mash,
        mash_k=args.mash_k, mash_s=args.mash_sketch_size,
        output_all_scores=args.output_all_scores, cache_dir=args.cache_dir,
        keep_intermediate=args.keep_intermediate, threads=args.threads,
        overwrite=args.overwrite, log_file=args.log_file, quiet=args.quiet,
    )
    return 0


def _cmd_annotate(args) -> int:
    run_annotation(
        input_fastas=args.input, out_dir=args.out_dir,
        with_amr=args.with_amr, with_replicon=args.with_replicon,
        with_mobility=args.with_mobility, amr_db=args.amr_db,
        abricate_minid=args.abricate_minid, abricate_mincov=args.abricate_mincov,
        prodigal_g=args.prodigal_g, prodigal_quiet=args.prodigal_quiet,
        macsy_evalue_search=args.macsy_evalue_search,
        macsy_i_evalue_sel=args.macsy_i_evalue_sel,
        macsy_coverage_profile=args.macsy_coverage_profile,
        macsy_no_cut_ga=args.macsy_no_cut_ga, macsy_force=args.macsy_force,
        models_dir=args.models_dir, threads=args.threads, overwrite=args.overwrite,
        log_file=args.log_file, quiet=args.quiet,
    )
    return 0


def _cmd_run(args) -> int:
    run_pipeline(
        input_fastas=args.input, out_dir=args.out_dir, model_dir=args.model_dir,
        database_dir=args.database, topk=args.topk, score_delta=args.score_delta,
        tax_level=args.tax_level, range_threshold_narrow=args.range_threshold_narrow,
        range_threshold_broad=args.range_threshold_broad,
        ptu_threshold=args.ptu_threshold, ptu_mash=args.ptu_mash,
        mash_k=args.mash_k, mash_s=args.mash_sketch_size,
        with_amr=args.with_amr, with_replicon=args.with_replicon,
        with_mobility=args.with_mobility, amr_db=args.amr_db,
        abricate_minid=args.abricate_minid, abricate_mincov=args.abricate_mincov,
        prodigal_g=args.prodigal_g, prodigal_quiet=args.prodigal_quiet,
        macsy_evalue_search=args.macsy_evalue_search,
        macsy_i_evalue_sel=args.macsy_i_evalue_sel,
        macsy_coverage_profile=args.macsy_coverage_profile,
        macsy_no_cut_ga=args.macsy_no_cut_ga, macsy_force=args.macsy_force,
        models_dir=args.models_dir, output_all_scores=args.output_all_scores,
        cache_dir=args.cache_dir, keep_intermediate=args.keep_intermediate,
        threads=args.threads, overwrite=args.overwrite, log_file=args.log_file,
        quiet=args.quiet,
    )
    return 0


def _cmd_version(args) -> int:
    _ = args
    print_header()
    print_info(f"PlasRiP {__version__}")
    print_info("Core modules: host-range prediction, PTU/PGG assignment")
    print_info("Reference architecture: external versioned database")
    return 0


def _database_report(database: str | None, deep: bool = False) -> dict[str, object]:
    return validate_database_dir(database, check_sketches=True, check_all_sketches=deep)


def _cmd_models(args) -> int:
    report: dict[str, object] = {
        "plasrip_version": __version__,
        "dependencies": check_dependencies(),
    }
    try:
        report["model_dir"] = validate_model_dir(args.model_dir)
    except Exception as exc:
        report["model_dir_error"] = str(exc)
    if args.database or __import__("os").environ.get("PLASRIP_DATABASE"):
        try:
            report["database"] = _database_report(args.database, deep=False)
        except Exception as exc:
            report["database_error"] = str(exc)

    if args.json:
        print(json.dumps(report, indent=2, ensure_ascii=False, default=_json_default))
    else:
        print_header()
        print_info(f"PlasRiP {__version__}")
        print_info("Dependencies:")
        for name, ok in report["dependencies"].items():
            print(f"  - {name}: {'OK' if ok else 'MISSING'}")
        if "model_dir" in report:
            print_info(f"Host-range model: {report['model_dir']['root']}")
        elif "model_dir_error" in report:
            print(f"Model directory error: {report['model_dir_error']}")
        if "database" in report:
            d = report["database"]
            print_info(f"Reference database: {d['root']}")
            print(f"  - version: {d['database_version']}")
            print(f"  - reference plasmids: {d['reference_plasmids']}")
            print(f"  - PTUs: {d['PTUs']}")
            print(f"  - PGGs: {d['PGGs']}")
        elif "database_error" in report:
            print(f"Reference database error: {report['database_error']}")
    return 0


def _cmd_database(args) -> int:
    # ``--check`` accepts an explicit database root or, when given without a
    # value, resolves PLASRIP_DATABASE through the database interface.
    database = None if args.check == "__PLASRIP_ENV__" else args.check
    report = _database_report(database, deep=args.deep)
    if args.json:
        print(json.dumps(report, indent=2, ensure_ascii=False, default=_json_default))
    else:
        print_header()
        print_info("Reference database validation: OK")
        for k, v in report.items():
            print(f"  - {k}: {v}")
    return 0


def _add_common_predict_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input plasmid FASTA files")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("--topk", type=int, default=10)
    parser.add_argument("--score-delta", type=float, default=DEFAULT_SCORE_DELTA)
    parser.add_argument("--tax-level", choices=["species","genus","family"], default="species")
    parser.add_argument("--range-threshold-narrow", type=float, default=None)
    parser.add_argument("--range-threshold-broad", type=float, default=None)
    parser.add_argument("--output-all-scores", action="store_true")
    parser.add_argument("--cache-dir", default=None)
    parser.add_argument("--keep-intermediate", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--quiet", action="store_true")


def _add_database_arg(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--database",
        metavar="DIR",
        default=None,
        help="PlasRiP reference-database root (or set PLASRIP_DATABASE)",
    )


def _add_database_ptu_args(parser: argparse.ArgumentParser) -> None:
    _add_database_arg(parser)
    parser.add_argument(
        "--ptu-threshold",
        type=float,
        default=None,
        help="PTU Mash-distance cutoff; default is read from database_manifest.tsv",
    )
    parser.add_argument("--ptu-mash", default="mash", help="Mash executable for PTU assignment")
    parser.add_argument("--mash-k", type=int, default=None, help="Advanced override; must match database manifest")
    parser.add_argument(
        "--mash-sketch-size",
        type=int,
        default=None,
        help="Advanced override; must match database manifest",
    )


def _add_annotation_args(parser: argparse.ArgumentParser) -> None:
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


def add_run_parser(subparsers) -> None:
    parser = subparsers.add_parser("run", help="Run full workflow")
    _add_common_predict_args(parser)
    parser.add_argument("--model-dir", default=None)
    _add_database_ptu_args(parser)
    _add_annotation_args(parser)
    parser.set_defaults(func=_cmd_run)


def add_predict_parser(subparsers) -> None:
    parser = subparsers.add_parser("predict", help="Predict plasmid host range")
    _add_common_predict_args(parser)
    parser.add_argument("--model-dir", default=None)
    parser.set_defaults(func=_cmd_predict)


def add_assign_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "assign",
        help="Assign plasmids to stable PTUs and report PGG context",
    )
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input plasmid FASTA files")
    _add_database_arg(parser)
    parser.add_argument("-o", "--output", required=True, help="PTU/PGG assignment TSV")
    parser.add_argument(
        "--ptu-threshold",
        type=float,
        default=None,
        help="PTU Mash-distance cutoff; default is read from database_manifest.tsv",
    )
    parser.add_argument("-t", "--threads", type=int, default=1, help="Mash thread count")
    parser.add_argument("--mash", default="mash", help="Mash executable")
    parser.add_argument("--mash-k", type=int, default=None, help="Advanced override; must match database manifest")
    parser.add_argument(
        "--mash-sketch-size",
        type=int,
        default=None,
        help="Advanced override; must match database manifest",
    )
    parser.add_argument("--quiet", action="store_true")
    parser.set_defaults(func=_cmd_assign)


def add_core_parser(subparsers) -> None:
    parser = subparsers.add_parser("core", help="Run host-range prediction + PTU/PGG assignment")
    _add_common_predict_args(parser)
    parser.add_argument("--model-dir", default=None)
    _add_database_ptu_args(parser)
    parser.set_defaults(func=_cmd_core)


def add_annotate_parser(subparsers) -> None:
    parser = subparsers.add_parser("annotate", help="Run one-stop plasmid annotation")
    parser.add_argument("-i", "--input", nargs="+", required=True)
    parser.add_argument("-o", "--out-dir", required=True)
    _add_annotation_args(parser)
    parser.add_argument("-t", "--threads", type=int, default=1)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--log-file", default=None)
    parser.add_argument("--quiet", action="store_true")
    parser.set_defaults(func=_cmd_annotate)


def add_version_parser(subparsers) -> None:
    parser = subparsers.add_parser("version", help="Show version information")
    parser.set_defaults(func=_cmd_version)


def add_models_parser(subparsers) -> None:
    parser = subparsers.add_parser("models", help="Check model, dependency, and database availability")
    parser.add_argument("--model-dir", default=None, help="Host-range model directory")
    _add_database_arg(parser)
    parser.add_argument("--json", action="store_true")
    parser.set_defaults(func=_cmd_models)


def add_database_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "database",
        help="Validate and inspect a PlasRiP reference database",
        description="Validate a schema-2.0 PlasRiP reference database.",
    )
    parser.add_argument(
        "--check",
        metavar="DIR",
        nargs="?",
        const="__PLASRIP_ENV__",
        required=True,
        help="Database root to validate; omit DIR to use PLASRIP_DATABASE",
    )
    parser.add_argument(
        "--deep",
        action="store_true",
        help="Validate every reference .msh file instead of representative sketches only",
    )
    parser.add_argument("--json", action="store_true", help="Emit the validation report as JSON")
    parser.set_defaults(func=_cmd_database)


class PlasRipArgumentParser(argparse.ArgumentParser):
    def format_help(self) -> str:
        return HEADER_ART + super().format_help()


def build_parser() -> argparse.ArgumentParser:
    parser = PlasRipArgumentParser(prog="plasrip", description="PlasRiP v0.2 command line interface")
    subparsers = parser.add_subparsers(dest="command")
    add_run_parser(subparsers)
    add_predict_parser(subparsers)
    add_assign_parser(subparsers)
    add_core_parser(subparsers)
    add_annotate_parser(subparsers)
    add_version_parser(subparsers)
    add_models_parser(subparsers)
    add_database_parser(subparsers)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    func = getattr(args, "func", None)
    if func is None:
        parser.print_help()
        return 0
    try:
        return int(func(args) or 0)
    except PlasRipError as exc:
        print_error(str(exc))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
