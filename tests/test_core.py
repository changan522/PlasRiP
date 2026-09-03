from __future__ import annotations

from types import SimpleNamespace

import pandas as pd

from plasrip.core import _ptu_annotation_table, build_core_summary
from plasrip.schemas import CORE_SUMMARY_COLUMNS


def _fake_db():
    ptu = pd.DataFrame([
        {
            "PTU_ID": "PTU000001",
            "representative_plasmid": "rep1",
            "n_plasmids": "38",
            "PGG_ID": "PGG-000005",
            "PGG_assignment_status": "ASSIGNED_STABLE_PGG",
            "mean_length_bp": "100000",
            "observed_host_range_level": "3",
            "observed_host_range_rank": "Family",
            "dominant_taxon": "Enterobacteriaceae",
            "dominant_taxon_coverage": "1",
            "n_host_species": "8",
            "top_host_species": "Enterobacter hormaechei",
            "top_host_species_coverage": "0.5",
            "n_host_genera": "4",
            "top_host_genus": "Enterobacter",
            "top_host_genus_coverage": "0.552632",
            "n_host_phyla": "1",
            "top_host_phylum": "Pseudomonadota",
            "top_host_phylum_coverage": "1",
            "dominant_replicon": "IncFIB",
            "dominant_replicon_fraction": "0.9",
            "dominant_mobility": "CONJ",
            "ARG_positive_fraction": "0.75",
            "n_unique_ARGs": "11",
            "unique_ARG_classes": "Beta-lactams;Aminoglycosides",
            "virulence_positive_fraction": "0.25",
            "n_unique_virulence_genes": "4",
        },
        {
            "PTU_ID": "PTU000003",
            "representative_plasmid": "rep3",
            "n_plasmids": "5",
            "PGG_ID": "",
            "PGG_assignment_status": "UNASSIGNED_OUTSIDE_CONSENSUS_SPACE",
            "mean_length_bp": "80000",
            "observed_host_range_level": "1",
            "observed_host_range_rank": "Species",
            "dominant_taxon": "Escherichia coli",
            "dominant_taxon_coverage": "1",
            "n_host_species": "1",
            "top_host_species": "Escherichia coli",
            "top_host_species_coverage": "1",
            "n_host_genera": "1",
            "top_host_genus": "Escherichia",
            "top_host_genus_coverage": "1",
            "n_host_phyla": "1",
            "top_host_phylum": "Pseudomonadota",
            "top_host_phylum_coverage": "1",
            "dominant_replicon": "IncFII",
            "dominant_replicon_fraction": "1",
            "dominant_mobility": "MOB",
            "ARG_positive_fraction": "0.2",
            "n_unique_ARGs": "1",
            "unique_ARG_classes": "Beta-lactams",
            "virulence_positive_fraction": "0",
            "n_unique_virulence_genes": "0",
        },
    ])
    pgg = pd.DataFrame([
        {
            "PGG_ID": "PGG-000005",
            "n_PTUs": "21",
            "n_plasmids": "300",
            "anchor_PTU": "PTU000051",
            "anchor_representative_plasmid": "rep51",
            "minimum_internal_support": "0.84375",
            "merge_with_parent_support": "0.5",
            "effective_persistence_width": "0.34375",
            "persistence_volume": "7.21875",
            "anchor_mean_within_group_Mash": "0.08",
        }
    ])
    return SimpleNamespace(ptu_metadata=ptu, pgg_metadata=pgg)


def _hostrange(plasmid="query1"):
    return pd.DataFrame([
        {
            "PLASMID": plasmid,
            "Predicted_host_top1": "Klebsiella pneumoniae",
            "Predicted_host_top1_score": 0.91,
            "Predicted_host_set": "Klebsiella pneumoniae;Escherichia coli",
            "Predicted_host_set_size": 2,
            "Predicted_host_set_cutoff_score": 0.71,
            "Predicted_host_score_delta": 0.20,
            "Predicted_genus_set": "Klebsiella;Escherichia",
            "Predicted_genus_set_size": 2,
            "Predicted_host_range_class": "broad",
        }
    ])


def test_step6_annotation_table_renames_observed_ptu_and_pgg_context():
    ann = _ptu_annotation_table(_fake_db())
    row = ann.loc[ann["PTU"] == "PTU000001"].iloc[0]

    assert row["PTU_observed_host_range_rank"] == "Family"
    assert row["PTU_observed_top_host_species"] == "Enterobacter hormaechei"
    assert row["PTU_observed_dominant_taxon"] == "Enterobacteriaceae"
    assert row["PTU_mean_length_bp"] == "100000"
    assert row["PTU_virulence_positive_fraction"] == "0.25"
    assert row["PGG_n_PTUs"] == "21"
    assert row["PGG_anchor_representative_plasmid"] == "rep51"
    assert row["PGG_effective_persistence_width"] == "0.34375"


def test_step6_core_summary_keeps_predicted_and_observed_host_information_distinct():
    ptu = pd.DataFrame([
        {
            "PLASMID": "query1",
            "PTU": "PTU000001",
            "PGG_ID": "PGG-000005",
            "PGG_assignment_status": "ASSIGNED_STABLE_PGG",
            "PTU_hits": 3,
            "how_assigned": "all_members",
            "closest_rep_dist": 0.012,
        }
    ])
    ann = _ptu_annotation_table(_fake_db())
    summary = build_core_summary(_hostrange(), ptu, ann)
    row = summary.iloc[0]

    assert list(summary.columns) == CORE_SUMMARY_COLUMNS
    assert row["Predicted_host_top1"] == "Klebsiella pneumoniae"
    assert row["Predicted_host_range_class"] == "broad"
    assert row["PTU_observed_host_range_rank"] == "Family"
    assert row["PTU_observed_top_host_species"] == "Enterobacter hormaechei"
    assert row["PTU_assignment_method"] == "all_members"
    assert row["PTU_candidate_hits"] == 3
    assert row["PTU_closest_rep_dist"] == 0.012
    assert row["PGG_ID"] == "PGG-000005"
    assert row["PGG_n_PTUs"] == "21"


def test_step6_ptu_without_stable_pgg_retains_reason_but_has_no_pgg_metrics():
    ptu = pd.DataFrame([
        {
            "PLASMID": "query2",
            "PTU": "PTU000003",
            "PGG_ID": "",
            "PGG_assignment_status": "UNASSIGNED_OUTSIDE_CONSENSUS_SPACE",
            "PTU_hits": 1,
            "how_assigned": "single_ptu_by_rep",
            "closest_rep_dist": 0.01,
        }
    ])
    ann = _ptu_annotation_table(_fake_db())
    summary = build_core_summary(_hostrange("query2"), ptu, ann)
    row = summary.iloc[0]

    assert row["PTU"] == "PTU000003"
    assert row["PGG_ID"] == ""
    assert row["PGG_assignment_status"] == "UNASSIGNED_OUTSIDE_CONSENSUS_SPACE"
    assert row["PTU_observed_host_range_rank"] == "Species"
    assert pd.isna(row["PGG_n_PTUs"])


def test_step6_no_ptu_does_not_inherit_reference_metadata():
    ptu = pd.DataFrame([
        {
            "PLASMID": "query3",
            "PTU": "NA",
            "PGG_ID": "NA",
            "PGG_assignment_status": "NOT_APPLICABLE_NO_PTU",
            "PTU_hits": 0,
            "how_assigned": "no_hit",
            "closest_rep_dist": None,
        }
    ])
    ann = _ptu_annotation_table(_fake_db())
    summary = build_core_summary(_hostrange("query3"), ptu, ann)
    row = summary.iloc[0]

    assert row["PTU"] == "NA"
    assert row["PGG_assignment_status"] == "NOT_APPLICABLE_NO_PTU"
    assert row["PTU_assignment_method"] == "no_hit"
    assert pd.isna(row["PTU_observed_host_range_rank"])
    assert pd.isna(row["PGG_n_PTUs"])


def test_step6_core_summary_can_be_built_when_host_prediction_row_is_missing():
    ptu = pd.DataFrame([
        {
            "PLASMID": "query4",
            "PTU": "PTU000001",
            "PGG_ID": "PGG-000005",
            "PGG_assignment_status": "ASSIGNED_STABLE_PGG",
            "PTU_hits": 1,
            "how_assigned": "closest_rep",
            "closest_rep_dist": 0.02,
        }
    ])
    ann = _ptu_annotation_table(_fake_db())
    summary = build_core_summary(pd.DataFrame(), ptu, ann)

    assert len(summary) == 1
    assert summary.iloc[0]["PLASMID"] == "query4"
    assert summary.iloc[0]["PTU"] == "PTU000001"
    assert pd.isna(summary.iloc[0]["Predicted_host_top1"])
