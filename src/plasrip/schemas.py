"""Column schemas and DataFrame helpers for PlasRiP outputs.

Public summary naming deliberately distinguishes three information layers:

1. ``Predicted_*``: host spectra predicted by the PlasRiP model for the query plasmid.
2. ``MPC_*`` / ``MPC_observed_*``: properties observed among reference plasmids in the
   assigned MPC. These are database annotations, not model predictions for the query.
3. ``PGG_*``: higher-order reference-group context for the assigned MPC.

The standalone ``mpc.tsv`` schema remains frozen by Step 5 and therefore retains its
compact assignment field names (for example ``how_assigned``). The combined
``core_summary.tsv`` uses explicit ``MPC_assignment_*`` names to avoid ambiguity.
"""

from __future__ import annotations

from typing import Iterable, Sequence
import pandas as pd

HOSTRANGE_COLUMNS = [
    "PLASMID",
    "Predicted_host_top1",
    "Predicted_host_top1_score",
    "Predicted_host_set",
    "Predicted_host_set_size",
    "Predicted_host_set_cutoff_score",
    "Predicted_host_score_delta",
    "Predicted_genus_set",
    "Predicted_genus_set_size",
    "Predicted_host_range_class",
]

HOSTRANGE_TOP_HITS_COLUMNS = [
    "PLASMID","Rank","Host_id","Host_species","Host_genus","Host_family","Score",
]

HOSTRANGE_SET_COLUMNS = [
    "PLASMID","Set_rank","Host_id","Host_species","Host_genus","Host_family","Score",
    "Predicted_host_top1_score","Predicted_host_set_cutoff_score","Predicted_host_score_delta",
]

HOSTRANGE_TAXONOMY_COLUMNS = [
    "PLASMID","Tax_level","Taxon","Aggregated_score","Rank","N_hosts",
]

# Step-5 frozen standalone MPC/PGG assignment output.
MPC_COLUMNS = [
    "PLASMID",
    "MPC",
    "PGG_ID",
    "PGG_assignment_status",
    "MPC_hits",
    "how_assigned",
    "closest_rep_dist",
    "MPC_closest",
    "MPC_closest_rep",
    "MPC_second",
    "MPC_second_rep_dist",
    "MPC_rep_dist_gap",
    "MPC_all_member",
    "MPC_all_member_max_dist",
    "MPC_all_member_avg_dist",
    "MPC_all_member_n",
    "MPC_network",
    "MPC_network_coverage",
    "MPC_network_mean_dist_leq_threshold",
    "MPC_network_n_within_threshold",
]

MPC_HITS_COLUMNS = [
    "PLASMID",
    "MPC_ID",
    "PGG_ID",
    "PGG_assignment_status",
    "Representative",
    "dist",
]

# Explicit assignment-context names used only in combined summaries.
MPC_ASSIGNMENT_SUMMARY_COLUMNS = [
    "MPC",
    "PGG_ID",
    "PGG_assignment_status",
    "MPC_assignment_method",
    "MPC_candidate_hits",
    "MPC_closest_rep_dist",
]

# Compact but information-complete MPC reference context for the combined summary.
# Full list-valued reference metadata remains in MPC_metadata.tsv in the database.
MPC_REFERENCE_METADATA_OUTPUT_COLUMNS = [
    "MPC_representative_plasmid",
    "MPC_n_plasmids",
    "MPC_mean_length_bp",
    "MPC_observed_host_range_level",
    "MPC_observed_host_range_rank",
    "MPC_observed_dominant_taxon",
    "MPC_observed_dominant_taxon_coverage",
    "MPC_observed_n_host_species",
    "MPC_observed_top_host_species",
    "MPC_observed_top_host_species_coverage",
    "MPC_observed_n_host_genera",
    "MPC_observed_top_host_genus",
    "MPC_observed_top_host_genus_coverage",
    "MPC_observed_n_host_phyla",
    "MPC_observed_top_host_phylum",
    "MPC_observed_top_host_phylum_coverage",
    "MPC_dominant_replicon",
    "MPC_dominant_replicon_fraction",
    "MPC_dominant_mobility",
    "MPC_ARG_positive_fraction",
    "MPC_n_unique_ARGs",
    "MPC_unique_ARG_classes",
    "MPC_virulence_positive_fraction",
    "MPC_n_unique_virulence_genes",
]

# PGG context is only populated when the assigned MPC belongs to a stable PGG.
PGG_REFERENCE_METADATA_OUTPUT_COLUMNS = [
    "PGG_n_MPCs",
    "PGG_n_plasmids",
    "PGG_anchor_MPC",
    "PGG_anchor_representative_plasmid",
    "PGG_minimum_internal_support",
    "PGG_merge_with_parent_support",
    "PGG_effective_persistence_width",
    "PGG_persistence_volume",
    "PGG_anchor_mean_within_group_Mash",
]

# Backward-compatible alias for internal imports created during early v0.2 development.
# New code should use MPC_REFERENCE_METADATA_OUTPUT_COLUMNS.
MPC_METADATA_OUTPUT_COLUMNS = MPC_REFERENCE_METADATA_OUTPUT_COLUMNS

CORE_SUMMARY_COLUMNS = (
    HOSTRANGE_COLUMNS
    + MPC_ASSIGNMENT_SUMMARY_COLUMNS
    + MPC_REFERENCE_METADATA_OUTPUT_COLUMNS
    + PGG_REFERENCE_METADATA_OUTPUT_COLUMNS
)

ANNOTATION_SUMMARY_COLUMNS = [
    "PLASMID","AMR_hits","AMR_genes","REPLICON_hits","REPLICON_genes","MOB",
]

# User-facing final report.  Keep this intentionally compact: detailed host-prediction
# diagnostics remain in hostrange.tsv, MPC/PGG reference context remains in
# core_summary.tsv, and assignment diagnostics remain in mpc.tsv.
SUMMARY_COLUMNS = [
    "PLASMID",

    # Query host prediction
    "Predicted_host_top1",
    "Predicted_host_top1_score",
    "Predicted_host_set",
    "Predicted_host_set_size",
    "Predicted_host_range_class",

    # MPC / PGG assignment
    "MPC",
    "PGG_ID",
    "PGG_assignment_status",
    "MPC_assignment_method",
    "MPC_candidate_hits",
    "MPC_closest_rep_dist",

    # Concise observed MPC context
    "MPC_n_plasmids",
    "MPC_observed_host_range_rank",
    "MPC_observed_dominant_taxon",
    "MPC_observed_n_host_species",
    "MPC_observed_top_host_species",
    "MPC_observed_top_host_species_coverage",
    "MPC_dominant_replicon",
    "MPC_dominant_mobility",

    # Query-plasmid annotation
    "AMR_hits",
    "AMR_genes",
    "REPLICON_hits",
    "REPLICON_genes",
    "MOB",
]


def ensure_columns(df: pd.DataFrame, columns: Sequence[str], fill_value: object = None) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        if col not in out.columns:
            out[col] = fill_value
    return out


def reorder_columns(df: pd.DataFrame, columns: Sequence[str], keep_extra: bool = True) -> pd.DataFrame:
    out = ensure_columns(df, columns)
    ordered = list(columns)
    if keep_extra:
        ordered.extend([c for c in out.columns if c not in ordered])
    return out.loc[:, ordered]


def empty_df(columns: Iterable[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=list(columns))


def empty_hostrange_df() -> pd.DataFrame:
    return empty_df(HOSTRANGE_COLUMNS)


def empty_hostrange_top_hits_df() -> pd.DataFrame:
    return empty_df(HOSTRANGE_TOP_HITS_COLUMNS)


def empty_hostrange_set_df() -> pd.DataFrame:
    return empty_df(HOSTRANGE_SET_COLUMNS)


def empty_hostrange_taxonomy_df() -> pd.DataFrame:
    return empty_df(HOSTRANGE_TAXONOMY_COLUMNS)


def empty_mpc_df() -> pd.DataFrame:
    return empty_df(MPC_COLUMNS)


def empty_mpc_hits_df() -> pd.DataFrame:
    return empty_df(MPC_HITS_COLUMNS)


def empty_core_summary_df() -> pd.DataFrame:
    return empty_df(CORE_SUMMARY_COLUMNS)


def empty_annotation_summary_df() -> pd.DataFrame:
    return empty_df(ANNOTATION_SUMMARY_COLUMNS)


def empty_summary_df() -> pd.DataFrame:
    return empty_df(SUMMARY_COLUMNS)
