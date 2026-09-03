"""Column schemas and DataFrame helpers for PlasRiP outputs.

Public summary naming deliberately distinguishes three information layers:

1. ``Predicted_*``: host spectra predicted by the PlasRiP model for the query plasmid.
2. ``PTU_*`` / ``PTU_observed_*``: properties observed among reference plasmids in the
   assigned PTU. These are database annotations, not model predictions for the query.
3. ``PGG_*``: higher-order reference-group context for the assigned PTU.

The standalone ``ptu.tsv`` schema remains frozen by Step 5 and therefore retains its
compact assignment field names (for example ``how_assigned``). The combined
``core_summary.tsv`` uses explicit ``PTU_assignment_*`` names to avoid ambiguity.
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

# Step-5 frozen standalone PTU/PGG assignment output.
PTU_COLUMNS = [
    "PLASMID",
    "PTU",
    "PGG_ID",
    "PGG_assignment_status",
    "PTU_hits",
    "how_assigned",
    "closest_rep_dist",
    "PTU_closest",
    "PTU_closest_rep",
    "PTU_second",
    "PTU_second_rep_dist",
    "PTU_rep_dist_gap",
    "PTU_all_member",
    "PTU_all_member_max_dist",
    "PTU_all_member_avg_dist",
    "PTU_all_member_n",
    "PTU_network",
    "PTU_network_coverage",
    "PTU_network_mean_dist_leq_threshold",
    "PTU_network_n_within_threshold",
]

PTU_HITS_COLUMNS = [
    "PLASMID",
    "PTU_ID",
    "PGG_ID",
    "PGG_assignment_status",
    "Representative",
    "dist",
]

# Explicit assignment-context names used only in combined summaries.
PTU_ASSIGNMENT_SUMMARY_COLUMNS = [
    "PTU",
    "PGG_ID",
    "PGG_assignment_status",
    "PTU_assignment_method",
    "PTU_candidate_hits",
    "PTU_closest_rep_dist",
]

# Compact but information-complete PTU reference context for the combined summary.
# Full list-valued reference metadata remains in PTU_metadata.tsv in the database.
PTU_REFERENCE_METADATA_OUTPUT_COLUMNS = [
    "PTU_representative_plasmid",
    "PTU_n_plasmids",
    "PTU_mean_length_bp",
    "PTU_observed_host_range_level",
    "PTU_observed_host_range_rank",
    "PTU_observed_dominant_taxon",
    "PTU_observed_dominant_taxon_coverage",
    "PTU_observed_n_host_species",
    "PTU_observed_top_host_species",
    "PTU_observed_top_host_species_coverage",
    "PTU_observed_n_host_genera",
    "PTU_observed_top_host_genus",
    "PTU_observed_top_host_genus_coverage",
    "PTU_observed_n_host_phyla",
    "PTU_observed_top_host_phylum",
    "PTU_observed_top_host_phylum_coverage",
    "PTU_dominant_replicon",
    "PTU_dominant_replicon_fraction",
    "PTU_dominant_mobility",
    "PTU_ARG_positive_fraction",
    "PTU_n_unique_ARGs",
    "PTU_unique_ARG_classes",
    "PTU_virulence_positive_fraction",
    "PTU_n_unique_virulence_genes",
]

# PGG context is only populated when the assigned PTU belongs to a stable PGG.
PGG_REFERENCE_METADATA_OUTPUT_COLUMNS = [
    "PGG_n_PTUs",
    "PGG_n_plasmids",
    "PGG_anchor_PTU",
    "PGG_anchor_representative_plasmid",
    "PGG_minimum_internal_support",
    "PGG_merge_with_parent_support",
    "PGG_effective_persistence_width",
    "PGG_persistence_volume",
    "PGG_anchor_mean_within_group_Mash",
]

# Backward-compatible alias for internal imports created during early v0.2 development.
# New code should use PTU_REFERENCE_METADATA_OUTPUT_COLUMNS.
PTU_METADATA_OUTPUT_COLUMNS = PTU_REFERENCE_METADATA_OUTPUT_COLUMNS

CORE_SUMMARY_COLUMNS = (
    HOSTRANGE_COLUMNS
    + PTU_ASSIGNMENT_SUMMARY_COLUMNS
    + PTU_REFERENCE_METADATA_OUTPUT_COLUMNS
    + PGG_REFERENCE_METADATA_OUTPUT_COLUMNS
)

ANNOTATION_SUMMARY_COLUMNS = [
    "PLASMID","AMR_hits","AMR_genes","REPLICON_hits","REPLICON_genes","MOB",
]

# User-facing final report.  Keep this intentionally compact: detailed host-prediction
# diagnostics remain in hostrange.tsv, PTU/PGG reference context remains in
# core_summary.tsv, and assignment diagnostics remain in ptu.tsv.
SUMMARY_COLUMNS = [
    "PLASMID",

    # Query host prediction
    "Predicted_host_top1",
    "Predicted_host_top1_score",
    "Predicted_host_set",
    "Predicted_host_set_size",
    "Predicted_host_range_class",

    # PTU / PGG assignment
    "PTU",
    "PGG_ID",
    "PGG_assignment_status",
    "PTU_assignment_method",
    "PTU_candidate_hits",
    "PTU_closest_rep_dist",

    # Concise observed PTU context
    "PTU_n_plasmids",
    "PTU_observed_host_range_rank",
    "PTU_observed_dominant_taxon",
    "PTU_observed_n_host_species",
    "PTU_observed_top_host_species",
    "PTU_observed_top_host_species_coverage",
    "PTU_dominant_replicon",
    "PTU_dominant_mobility",

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


def empty_ptu_df() -> pd.DataFrame:
    return empty_df(PTU_COLUMNS)


def empty_ptu_hits_df() -> pd.DataFrame:
    return empty_df(PTU_HITS_COLUMNS)


def empty_core_summary_df() -> pd.DataFrame:
    return empty_df(CORE_SUMMARY_COLUMNS)


def empty_annotation_summary_df() -> pd.DataFrame:
    return empty_df(ANNOTATION_SUMMARY_COLUMNS)


def empty_summary_df() -> pd.DataFrame:
    return empty_df(SUMMARY_COLUMNS)
