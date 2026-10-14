from __future__ import annotations

import pandas as pd

from plasrip import schemas


def test_empty_hostrange_df_columns():
    df = schemas.empty_hostrange_df()
    assert list(df.columns) == schemas.HOSTRANGE_COLUMNS
    assert df.empty


def test_reorder_columns_keeps_extra_columns():
    df = pd.DataFrame([{"extra": 1, "PLASMID": "p1"}])
    out = schemas.reorder_columns(df, schemas.SUMMARY_COLUMNS)
    assert list(out.columns[: len(schemas.SUMMARY_COLUMNS)]) == schemas.SUMMARY_COLUMNS
    assert out.loc[0, "extra"] == 1


def test_step5_mpc_columns_prioritize_mpc_and_pgg_identifiers():
    assert schemas.MPC_COLUMNS[:6] == [
        "PLASMID",
        "MPC",
        "PGG_ID",
        "PGG_assignment_status",
        "MPC_hits",
        "how_assigned",
    ]


def test_step5_mpc_hits_schema_includes_pgg_status():
    assert schemas.MPC_HITS_COLUMNS == [
        "PLASMID",
        "MPC_ID",
        "PGG_ID",
        "PGG_assignment_status",
        "Representative",
        "dist",
    ]


def test_step6_core_summary_layers_are_explicitly_separated():
    predicted = [c for c in schemas.CORE_SUMMARY_COLUMNS if c.startswith("Predicted_")]
    observed = [c for c in schemas.CORE_SUMMARY_COLUMNS if c.startswith("MPC_observed_")]
    pgg = [c for c in schemas.CORE_SUMMARY_COLUMNS if c.startswith("PGG_")]

    assert predicted
    assert observed
    assert pgg
    assert "MPC_assignment_method" in schemas.CORE_SUMMARY_COLUMNS
    assert "MPC_candidate_hits" in schemas.CORE_SUMMARY_COLUMNS
    assert "MPC_closest_rep_dist" in schemas.CORE_SUMMARY_COLUMNS
    assert "how_assigned" not in schemas.CORE_SUMMARY_COLUMNS
    assert "MPC_hits" not in schemas.CORE_SUMMARY_COLUMNS
    assert "closest_rep_dist" not in schemas.CORE_SUMMARY_COLUMNS


def test_step6_standalone_mpc_schema_remains_step5_frozen():
    assert "how_assigned" in schemas.MPC_COLUMNS
    assert "MPC_hits" in schemas.MPC_COLUMNS
    assert "closest_rep_dist" in schemas.MPC_COLUMNS
    assert "MPC_assignment_method" not in schemas.MPC_COLUMNS
