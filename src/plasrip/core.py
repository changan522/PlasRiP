"""Host-range prediction plus database-backed PTU/PGG assignment."""

from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd

from .assign import run_ptu_assignment
from .config import (
    DEFAULT_CORE_LOG_NAME,
    DEFAULT_CORE_SUMMARY_OUT,
    DEFAULT_HOSTRANGE_OUT,
    DEFAULT_HOSTRANGE_SET_OUT,
    DEFAULT_HOSTRANGE_TAXONOMY_OUT,
    DEFAULT_HOSTRANGE_TOP_HITS_OUT,
    DEFAULT_PTU_HITS_SUFFIX,
    DEFAULT_PTU_OUT,
)
from .database import build_assignment_indexes, load_database, ptu_metadata_lookup
from .predict import run_host_prediction
from .schemas import (
    CORE_SUMMARY_COLUMNS,
    PTU_COLUMNS,
    PTU_HITS_COLUMNS,
    PTU_ASSIGNMENT_SUMMARY_COLUMNS,
    PTU_REFERENCE_METADATA_OUTPUT_COLUMNS,
    PGG_REFERENCE_METADATA_OUTPUT_COLUMNS,
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


def _standardize_ptu_df(ptu_df: pd.DataFrame) -> pd.DataFrame:
    if ptu_df is None or ptu_df.empty:
        return empty_ptu_df()
    out = ptu_df.copy().rename(columns={
        "query": "PLASMID",
        "n_ptu_hits": "PTU_hits",
        "ptu_assigned": "PTU",
        "closest_ptu": "PTU_closest",
        "closest_rep": "PTU_closest_rep",
        "second_ptu": "PTU_second",
        "second_rep_dist": "PTU_second_rep_dist",
        "rep_dist_gap": "PTU_rep_dist_gap",
        "all_member_ptu": "PTU_all_member",
        "all_member_max_dist": "PTU_all_member_max_dist",
        "all_member_avg_dist": "PTU_all_member_avg_dist",
        "all_member_n": "PTU_all_member_n",
        "network_ptu": "PTU_network",
        "network_coverage": "PTU_network_coverage",
        "network_mean_dist_leq_threshold": "PTU_network_mean_dist_leq_threshold",
        "network_n_within_threshold": "PTU_network_n_within_threshold",
    })
    if "PLASMID" in out.columns:
        out["PLASMID"] = out["PLASMID"].map(_norm_plasmid_key)
    return reorder_columns(out, PTU_COLUMNS, keep_extra=False)


def _standardize_ptu_hits_df(ptu_hits_df: pd.DataFrame) -> pd.DataFrame:
    if ptu_hits_df is None or ptu_hits_df.empty:
        return empty_ptu_hits_df()
    out = ptu_hits_df.copy().rename(columns={"plasmid": "PLASMID"})
    if "PLASMID" in out.columns:
        out["PLASMID"] = out["PLASMID"].map(_norm_plasmid_key)
    return reorder_columns(out, PTU_HITS_COLUMNS, keep_extra=False)


def _ptu_annotation_table(db) -> pd.DataFrame:
    """Build the reference-context table merged into ``core_summary.tsv``.

    Naming is intentionally explicit:
    - ``Predicted_*`` fields elsewhere in the core summary describe model predictions
      for the query plasmid.
    - ``PTU_observed_*`` fields here describe host distributions observed among
      reference plasmids belonging to the assigned PTU.
    - ``PGG_*`` fields describe the higher-order reference group of that PTU.

    The database's list-valued host/ARG/virulence fields remain available in the
    external database and are deliberately not duplicated into the compact summary.
    """
    ptu = db.ptu_metadata.copy()
    ptu_rename = {
        "PTU_ID": "PTU",
        "representative_plasmid": "PTU_representative_plasmid",
        "n_plasmids": "PTU_n_plasmids",
        "mean_length_bp": "PTU_mean_length_bp",
        "observed_host_range_level": "PTU_observed_host_range_level",
        "observed_host_range_rank": "PTU_observed_host_range_rank",
        "dominant_taxon": "PTU_observed_dominant_taxon",
        "dominant_taxon_coverage": "PTU_observed_dominant_taxon_coverage",
        "n_host_species": "PTU_observed_n_host_species",
        "top_host_species": "PTU_observed_top_host_species",
        "top_host_species_coverage": "PTU_observed_top_host_species_coverage",
        "n_host_genera": "PTU_observed_n_host_genera",
        "top_host_genus": "PTU_observed_top_host_genus",
        "top_host_genus_coverage": "PTU_observed_top_host_genus_coverage",
        "n_host_phyla": "PTU_observed_n_host_phyla",
        "top_host_phylum": "PTU_observed_top_host_phylum",
        "top_host_phylum_coverage": "PTU_observed_top_host_phylum_coverage",
        "dominant_replicon": "PTU_dominant_replicon",
        "dominant_replicon_fraction": "PTU_dominant_replicon_fraction",
        "dominant_mobility": "PTU_dominant_mobility",
        "ARG_positive_fraction": "PTU_ARG_positive_fraction",
        "n_unique_ARGs": "PTU_n_unique_ARGs",
        "unique_ARG_classes": "PTU_unique_ARG_classes",
        "virulence_positive_fraction": "PTU_virulence_positive_fraction",
        "n_unique_virulence_genes": "PTU_n_unique_virulence_genes",
    }
    ptu_keep = [c for c in ptu_rename if c in ptu.columns] + [
        c for c in ["PGG_ID", "PGG_assignment_status"] if c in ptu.columns
    ]
    ptu = ptu.loc[:, ptu_keep].rename(columns=ptu_rename)

    pgg = db.pgg_metadata.copy()
    pgg_rename = {
        "n_PTUs": "PGG_n_PTUs",
        "n_plasmids": "PGG_n_plasmids",
        "anchor_PTU": "PGG_anchor_PTU",
        "anchor_representative_plasmid": "PGG_anchor_representative_plasmid",
        "minimum_internal_support": "PGG_minimum_internal_support",
        "merge_with_parent_support": "PGG_merge_with_parent_support",
        "effective_persistence_width": "PGG_effective_persistence_width",
        "persistence_volume": "PGG_persistence_volume",
        "anchor_mean_within_group_Mash": "PGG_anchor_mean_within_group_Mash",
    }
    pgg_keep = ["PGG_ID"] + [c for c in pgg_rename if c in pgg.columns]
    pgg = pgg.loc[:, pgg_keep].rename(columns=pgg_rename)

    annotation = ptu.merge(pgg, on="PGG_ID", how="left")
    expected = [
        "PTU", "PGG_ID", "PGG_assignment_status",
        *PTU_REFERENCE_METADATA_OUTPUT_COLUMNS,
        *PGG_REFERENCE_METADATA_OUTPUT_COLUMNS,
    ]
    return reorder_columns(annotation, expected, keep_extra=False)


def build_core_summary(
    hostrange_df: pd.DataFrame,
    ptu_df: pd.DataFrame,
    ptu_annotation: pd.DataFrame,
) -> pd.DataFrame:
    """Merge query predictions, PTU assignment, and PTU/PGG reference context.

    ``hostrange_df`` is the query-level model output and therefore retains the
    ``Predicted_*`` prefix. PTU assignment fields are renamed to explicit
    ``PTU_assignment_*`` names in the combined summary, while the standalone
    ``ptu.tsv`` remains unchanged according to the Step-5 schema.
    """
    predicted = hostrange_df.copy() if hostrange_df is not None else pd.DataFrame()
    assignment = ptu_df.copy() if ptu_df is not None else pd.DataFrame()

    if "PLASMID" in predicted.columns:
        predicted["PLASMID"] = predicted["PLASMID"].map(_norm_plasmid_key)
    if "PLASMID" in assignment.columns:
        assignment["PLASMID"] = assignment["PLASMID"].map(_norm_plasmid_key)

    if predicted.empty or "PLASMID" not in predicted.columns:
        if not assignment.empty and "PLASMID" in assignment.columns:
            predicted = pd.DataFrame({
                "PLASMID": assignment["PLASMID"].dropna().astype(str).drop_duplicates()
            })
        else:
            predicted = pd.DataFrame(columns=["PLASMID"])

    assignment_rename = {
        "how_assigned": "PTU_assignment_method",
        "PTU_hits": "PTU_candidate_hits",
        "closest_rep_dist": "PTU_closest_rep_dist",
    }
    assignment_cols = [
        "PLASMID", "PTU", "PGG_ID", "PGG_assignment_status",
        "how_assigned", "PTU_hits", "closest_rep_dist",
    ]
    assignment_sub = assignment.loc[:, [c for c in assignment_cols if c in assignment.columns]].copy()
    assignment_sub.rename(columns=assignment_rename, inplace=True)

    summary = predicted.copy()
    if not assignment_sub.empty and "PLASMID" in assignment_sub.columns:
        summary = summary.merge(assignment_sub, on="PLASMID", how="left")

    if not ptu_annotation.empty and "PTU" in summary.columns:
        # PTU assignment is authoritative for PTU/PGG identifiers and PGG status.
        reference_context = ptu_annotation.drop(
            columns=["PGG_ID", "PGG_assignment_status"], errors="ignore"
        )
        summary = summary.merge(reference_context, on="PTU", how="left")

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
    database_dir: str | None,
    out_dir: str,
    topk: int = 10,
    score_delta: float = 0.20,
    tax_level: str = "species",
    range_threshold_narrow: float | None = None,
    range_threshold_broad: float | None = None,
    ptu_threshold: float | None = None,
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
    _database=None,
    _assignment_indexes=None,
    _ptu_meta_lookup=None,
    _ptu_annotation=None,
) -> dict:
    out_root = ensure_dir(out_dir, overwrite=overwrite)
    logger = setup_logger("plasrip.core", log_file or str(out_root / DEFAULT_CORE_LOG_NAME), quiet=True)

    # Public calls load the external database once for the complete core run.
    # The full pipeline may provide the already-loaded database and derived
    # indexes through private keyword arguments so that multiple query plasmids
    # share one validated in-memory database object.
    db = _database if _database is not None else load_database(database_dir, check_sketches=True)
    assignment_indexes = (
        _assignment_indexes if _assignment_indexes is not None else build_assignment_indexes(db)
    )
    ptu_lookup = _ptu_meta_lookup if _ptu_meta_lookup is not None else ptu_metadata_lookup(db)
    ptu_annotation = _ptu_annotation if _ptu_annotation is not None else _ptu_annotation_table(db)

    frames = {k: [] for k in ["hostrange","top_hits","taxonomy","set_df","ptu","ptu_hits","core_summary","all_scores"]}
    total = len(input_fastas)
    completed = failed = 0

    try:
        _reset_core_outputs(out_root, output_all_scores=output_all_scores)
        if not quiet:
            if show_header:
                print_header()
            print_info(f"reference database : {db.database_version}")
            print_info(f"PTUs / PGGs        : {len(db.ptu_metadata)} / {len(db.pgg_metadata)}")

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
                    input_fastas=[fasta], model_dir=model_dir,
                    out_dir=str(workdir / "predict"), topk=topk, score_delta=score_delta,
                    tax_level=tax_level, range_threshold_narrow=range_threshold_narrow,
                    range_threshold_broad=range_threshold_broad,
                    output_all_scores=output_all_scores,
                    cache_dir=str(workdir / "cache" / "hostrange") if cache_dir else None,
                    keep_intermediate=keep_intermediate, threads=threads, overwrite=True,
                    log_file=None, quiet=True,
                )
                hr = host_results.get("hostrange", pd.DataFrame())
                if not hr.empty and not quiet:
                    row = hr.iloc[0]
                    print_done(
                        f"{name} :: top1={row.get('Predicted_host_top1', 'NA')} "
                        f"set_size={row.get('Predicted_host_set_size', 'NA')} "
                        f"range={row.get('Predicted_host_range_class', 'NA')}"
                    )

                if not quiet:
                    print_step(name, "PTU/PGG assignment")
                ptu_out = workdir / "ptu.tsv"
                ptu_raw = run_ptu_assignment(
                    input_fastas=[fasta], database_dir=str(db.root), output=str(ptu_out),
                    threshold=ptu_threshold, threads=threads, mash_cmd=ptu_mash,
                    mash_k=mash_k, mash_s=mash_s, show_progress=False,
                    _database=db, _assignment_indexes=assignment_indexes,
                    _ptu_meta_lookup=ptu_lookup,
                )
                ptu_df = _standardize_ptu_df(ptu_raw)
                hits_path = Path(str(ptu_out) + DEFAULT_PTU_HITS_SUFFIX)
                ptu_hits = (
                    _standardize_ptu_hits_df(pd.read_csv(hits_path, sep="\t", keep_default_na=False))
                    if hits_path.is_file() else empty_ptu_hits_df()
                )

                if not ptu_df.empty and not quiet:
                    prow = ptu_df.iloc[0]
                    print_done(
                        f"{name} :: PTU={prow.get('PTU', 'NA')} "
                        f"PGG={prow.get('PGG_ID') or 'NA'} "
                        f"method={prow.get('how_assigned', 'NA')}"
                    )

                results = dict(host_results)
                results["ptu"] = ptu_df
                results["ptu_hits"] = ptu_hits
                results["core_summary"] = build_core_summary(
                    results.get("hostrange", pd.DataFrame()), ptu_df, ptu_annotation
                )
                _append_core_outputs(results, out_root, output_all_scores=output_all_scores)

                for key in ["hostrange","top_hits","taxonomy","set_df","ptu","ptu_hits","core_summary"]:
                    frames[key].append(results.get(key, pd.DataFrame()))
                if output_all_scores and "all_scores" in results:
                    frames["all_scores"].append(results["all_scores"])
                completed += 1
            except Exception as exc:  # pragma: no cover - external tools/runtime
                failed += 1
                logger.exception("Failed plasmid: %s", fasta)
                if not quiet:
                    print_error(f"{name} :: {exc}")

        if not quiet:
            print_final_summary(total=total, completed=completed, failed=failed, out_dir=out_root)

        return {
            "hostrange": pd.concat(frames["hostrange"], ignore_index=True) if frames["hostrange"] else pd.DataFrame(),
            "top_hits": pd.concat(frames["top_hits"], ignore_index=True) if frames["top_hits"] else pd.DataFrame(),
            "taxonomy": pd.concat(frames["taxonomy"], ignore_index=True) if frames["taxonomy"] else pd.DataFrame(),
            "set_df": pd.concat(frames["set_df"], ignore_index=True) if frames["set_df"] else pd.DataFrame(),
            "ptu": pd.concat(frames["ptu"], ignore_index=True) if frames["ptu"] else empty_ptu_df(),
            "ptu_hits": pd.concat(frames["ptu_hits"], ignore_index=True) if frames["ptu_hits"] else empty_ptu_hits_df(),
            "core_summary": pd.concat(frames["core_summary"], ignore_index=True) if frames["core_summary"] else pd.DataFrame(columns=CORE_SUMMARY_COLUMNS),
            "all_scores": pd.concat(frames["all_scores"], ignore_index=True) if frames["all_scores"] else pd.DataFrame(),
        }
    finally:
        close_logger(logger)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run host-range prediction + PTU/PGG assignment")
    parser.add_argument("-i", "--input", nargs="+", required=True)
    parser.add_argument("--model-dir", default=None)
    parser.add_argument("--database", default=None, help="Reference database directory (or set PLASRIP_DATABASE)")
    parser.add_argument("-o", "--out-dir", required=True)
    parser.add_argument("--topk", type=int, default=10)
    parser.add_argument("--score-delta", type=float, default=0.20)
    parser.add_argument("--tax-level", choices=["species","genus","family"], default="species")
    parser.add_argument("--range-threshold-narrow", type=float, default=None)
    parser.add_argument("--range-threshold-broad", type=float, default=None)
    parser.add_argument("--ptu-threshold", type=float, default=None, help="Default: database manifest cutoff")
    parser.add_argument("--ptu-mash", default="mash")
    parser.add_argument("--mash-k", type=int, default=None, help="Must match database manifest")
    parser.add_argument("--mash-sketch-size", dest="mash_s", type=int, default=None, help="Must match database manifest")
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
        input_fastas=args.input, model_dir=args.model_dir, database_dir=args.database,
        out_dir=args.out_dir, topk=args.topk, score_delta=args.score_delta,
        tax_level=args.tax_level, range_threshold_narrow=args.range_threshold_narrow,
        range_threshold_broad=args.range_threshold_broad, ptu_threshold=args.ptu_threshold,
        ptu_mash=args.ptu_mash, mash_k=args.mash_k, mash_s=args.mash_s,
        output_all_scores=args.output_all_scores, cache_dir=args.cache_dir,
        keep_intermediate=args.keep_intermediate, threads=args.threads,
        overwrite=args.overwrite, log_file=args.log_file, quiet=args.quiet,
        show_header=True,
    )


if __name__ == "__main__":
    main()
