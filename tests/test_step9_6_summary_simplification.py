from __future__ import annotations

import pandas as pd

from plasrip import pipeline
from plasrip.schemas import CORE_SUMMARY_COLUMNS, SUMMARY_COLUMNS


EXPECTED_SUMMARY_COLUMNS = [
    "PLASMID",
    "Predicted_host_top1",
    "Predicted_host_top1_score",
    "Predicted_host_set",
    "Predicted_host_set_size",
    "Predicted_host_range_class",
    "MPC",
    "PGG_ID",
    "PGG_assignment_status",
    "MPC_assignment_method",
    "MPC_candidate_hits",
    "MPC_closest_rep_dist",
    "MPC_n_plasmids",
    "MPC_observed_host_range_rank",
    "MPC_observed_dominant_taxon",
    "MPC_observed_n_host_species",
    "MPC_observed_top_host_species",
    "MPC_observed_top_host_species_coverage",
    "MPC_dominant_replicon",
    "MPC_dominant_mobility",
    "AMR_hits",
    "AMR_genes",
    "REPLICON_hits",
    "REPLICON_genes",
    "MOB",
]


def test_final_summary_schema_is_compact_and_frozen():
    assert SUMMARY_COLUMNS == EXPECTED_SUMMARY_COLUMNS
    assert len(SUMMARY_COLUMNS) == 25


def test_detailed_fields_remain_in_core_summary_not_final_summary():
    detailed = {
        "Predicted_host_set_cutoff_score",
        "Predicted_host_score_delta",
        "Predicted_genus_set",
        "Predicted_genus_set_size",
        "MPC_representative_plasmid",
        "MPC_mean_length_bp",
        "MPC_observed_n_host_genera",
        "MPC_observed_top_host_genus",
        "MPC_observed_n_host_phyla",
        "MPC_observed_top_host_phylum",
        "MPC_ARG_positive_fraction",
        "MPC_virulence_positive_fraction",
        "PGG_n_MPCs",
        "PGG_anchor_MPC",
        "PGG_minimum_internal_support",
    }
    assert detailed <= set(CORE_SUMMARY_COLUMNS)
    assert detailed.isdisjoint(SUMMARY_COLUMNS)


def test_merge_summary_keeps_only_public_compact_schema():
    core = pd.DataFrame([{
        "PLASMID": "q1",
        "Predicted_host_top1": "Escherichia coli",
        "Predicted_host_top1_score": 0.91,
        "Predicted_host_set": "Escherichia coli;Klebsiella pneumoniae",
        "Predicted_host_set_size": 2,
        "Predicted_host_range_class": "intermediate",
        "Predicted_host_score_delta": 0.20,
        "MPC": "MPC000001",
        "PGG_ID": "PGG-000005",
        "PGG_assignment_status": "ASSIGNED_STABLE_PGG",
        "MPC_assignment_method": "closest_rep",
        "MPC_candidate_hits": 1,
        "MPC_closest_rep_dist": 0.01,
        "MPC_n_plasmids": 38,
        "MPC_observed_host_range_rank": "Order",
        "MPC_observed_dominant_taxon": "Enterobacterales",
        "MPC_observed_n_host_species": 8,
        "MPC_observed_top_host_species": "Enterobacter hormaechei",
        "MPC_observed_top_host_species_coverage": 0.5,
        "MPC_dominant_replicon": "IncFIB",
        "MPC_dominant_mobility": "CONJ",
        "PGG_n_MPCs": 21,
    }])
    ann = pd.DataFrame([{
        "PLASMID": "q1",
        "AMR_hits": 2,
        "AMR_genes": "blaNDM-5;blaTEM-1",
        "REPLICON_hits": 1,
        "REPLICON_genes": "IncFIB",
        "MOB": "CONJ",
    }])
    out = pipeline.build_final_summary(core, ann)
    assert list(out.columns) == EXPECTED_SUMMARY_COLUMNS
    assert out.loc[0, "MPC"] == "MPC000001"
    assert "Predicted_host_score_delta" not in out.columns
    assert "PGG_n_MPCs" not in out.columns
