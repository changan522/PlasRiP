"""Column schemas and DataFrame helpers for PlasRiP outputs."""

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
    "Predicted_host_confidence",
]

HOSTRANGE_TOP_HITS_COLUMNS = [
    "PLASMID",
    "Rank",
    "Host_id",
    "Host_species",
    "Host_genus",
    "Host_family",
    "Score",
]

HOSTRANGE_SET_COLUMNS = [
    "PLASMID",
    "Set_rank",
    "Host_id",
    "Host_species",
    "Host_genus",
    "Host_family",
    "Score",
    "Predicted_host_top1_score",
    "Predicted_host_set_cutoff_score",
    "Predicted_host_score_delta",
]

HOSTRANGE_TAXONOMY_COLUMNS = [
    "PLASMID",
    "Tax_level",
    "Taxon",
    "Aggregated_score",
    "Rank",
    "N_hosts",
]

PTU_COLUMNS = [
    "PLASMID",
    "PTU_hits",
    "PTU",
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
    "Representative",
    "dist",
]

PTU_INFO_COLUMNS = [
    "PTU",
    "PTU_n_plasmids",
    "PTU_host_range_level",
    "PTU_host_species_range",
    "PTU_mob_dominant",
    "PTU_avg_genome_size_bp",
    "PTU_genome_size_min_bp",
    "PTU_genome_size_max_bp",
    "PTU_genome_size_range_bp",
    "PTU_avg_AMR_genes_per_plasmid",
    "PTU_unique_AMR_classes",
]

CORE_SUMMARY_COLUMNS = [
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
    "Predicted_host_confidence",
    "PTU",
    "how_assigned",
    "PTU_hits",
    "closest_rep_dist",
    "PTU_n_plasmids",
    "PTU_host_range_level",
    "PTU_host_species_range",
    "PTU_mob_dominant",
    "PTU_avg_genome_size_bp",
    "PTU_genome_size_min_bp",
    "PTU_genome_size_max_bp",
    "PTU_genome_size_range_bp",
    "PTU_avg_AMR_genes_per_plasmid",
    "PTU_unique_AMR_classes",
]

ANNOTATION_SUMMARY_COLUMNS = [
    "PLASMID",
    "AMR_hits",
    "AMR_genes",
    "REPLICON_hits",
    "REPLICON_genes",
    "MOB",
]

SUMMARY_COLUMNS = [
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
    "Predicted_host_confidence",
    "PTU",
    "MOB",
    "how_assigned",
    "PTU_hits",
    "closest_rep_dist",
    "AMR_hits",
    "AMR_genes",
    "REPLICON_hits",
    "REPLICON_genes",
    "PTU_n_plasmids",
    "PTU_host_range_level",
    "PTU_host_species_range",
    "PTU_mob_dominant",
    "PTU_avg_genome_size_bp",
    "PTU_genome_size_min_bp",
    "PTU_genome_size_max_bp",
    "PTU_genome_size_range_bp",
    "PTU_avg_AMR_genes_per_plasmid",
    "PTU_unique_AMR_classes",
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
