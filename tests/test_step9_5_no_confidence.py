import pandas as pd

from plasrip.hostrange_postprocess import build_hostrange_summary
from plasrip.schemas import HOSTRANGE_COLUMNS, CORE_SUMMARY_COLUMNS, SUMMARY_COLUMNS


def test_prediction_confidence_removed_from_public_schemas():
    assert "Predicted_host_confidence" not in HOSTRANGE_COLUMNS
    assert "Predicted_host_confidence" not in CORE_SUMMARY_COLUMNS
    assert "Predicted_host_confidence" not in SUMMARY_COLUMNS


def test_hostrange_summary_reports_top1_score_without_confidence_label():
    top_hits = pd.DataFrame([
        {
            "PLASMID": "q1", "Rank": 1, "Host_id": "h1",
            "Host_species": "Species A", "Host_genus": "GenusA",
            "Host_family": "FamilyA", "Score": 0.87,
        },
        {
            "PLASMID": "q1", "Rank": 2, "Host_id": "h2",
            "Host_species": "Species B", "Host_genus": "GenusB",
            "Host_family": "FamilyB", "Score": 0.84,
        },
    ])
    set_df = pd.DataFrame([
        {
            "PLASMID": "q1", "Set_rank": 1, "Host_id": "h1",
            "Host_species": "Species A", "Host_genus": "GenusA",
            "Host_family": "FamilyA", "Score": 0.87,
            "Predicted_host_top1_score": 0.87,
            "Predicted_host_set_cutoff_score": 0.67,
            "Predicted_host_score_delta": 0.20,
        },
        {
            "PLASMID": "q1", "Set_rank": 2, "Host_id": "h2",
            "Host_species": "Species B", "Host_genus": "GenusB",
            "Host_family": "FamilyB", "Score": 0.84,
            "Predicted_host_top1_score": 0.87,
            "Predicted_host_set_cutoff_score": 0.67,
            "Predicted_host_score_delta": 0.20,
        },
    ])
    range_df = pd.DataFrame([
        {"PLASMID": "q1", "Predicted_host_range_class": "intermediate"}
    ])

    out = build_hostrange_summary(top_hits, set_df, range_df)
    assert list(out.columns) == HOSTRANGE_COLUMNS
    assert "Predicted_host_confidence" not in out.columns
    assert float(out.loc[0, "Predicted_host_top1_score"]) == 0.87
