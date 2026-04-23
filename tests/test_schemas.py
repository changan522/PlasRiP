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
