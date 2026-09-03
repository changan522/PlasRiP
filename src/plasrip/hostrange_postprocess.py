"""Postprocessing for PlasRiP host-range outputs aligned to the 3.2 training script."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

try:
    from .config import DEFAULT_SCORE_DELTA, DEFAULT_TOPK  # type: ignore
except Exception:  # pragma: no cover
    DEFAULT_TOPK = 10
    DEFAULT_SCORE_DELTA = 0.20

from .schemas import (
    HOSTRANGE_COLUMNS,
    HOSTRANGE_SET_COLUMNS,
    empty_hostrange_df,
    empty_hostrange_set_df,
    empty_hostrange_taxonomy_df,
    empty_hostrange_top_hits_df,
    reorder_columns,
)
from .utils.taxonomy import aggregate_scores_by_taxonomy, load_host_metadata


def _norm_plasmid_key(x: object) -> str:
    return Path(str(x)).name if x is not None else ""


def _first_unique_join(values: list[object]) -> str:
    out: list[str] = []
    seen: set[str] = set()
    for v in values:
        s = str(v or "").strip()
        if not s or s in seen:
            continue
        out.append(s)
        seen.add(s)
    return ";".join(out)


def build_top_hits(score_df: pd.DataFrame, host_meta_df: pd.DataFrame, topk: int = DEFAULT_TOPK) -> pd.DataFrame:
    if score_df is None or score_df.empty:
        return empty_hostrange_top_hits_df()
    out = score_df.copy()
    out["PLASMID"] = out["PLASMID"].map(_norm_plasmid_key)
    out = out.merge(host_meta_df, on="Host_id", how="left")
    out = out.sort_values(["PLASMID", "Score"], ascending=[True, False]).copy()
    out["Rank"] = out.groupby("PLASMID").cumcount() + 1
    out = out.loc[out["Rank"] <= int(max(1, topk))].copy()
    cols = ["PLASMID", "Rank", "Host_id", "Host_species", "Host_genus", "Host_family", "Score"]
    return reorder_columns(out.loc[:, cols], cols, keep_extra=False)


def build_prediction_sets(
    score_df: pd.DataFrame,
    host_meta_df: pd.DataFrame,
    score_delta: float = DEFAULT_SCORE_DELTA,
) -> pd.DataFrame:
    """Build the 3.2 prediction set using score >= top1_score - score_delta."""
    if score_df is None or score_df.empty:
        return empty_hostrange_set_df()

    out = score_df.copy()
    out["PLASMID"] = out["PLASMID"].map(_norm_plasmid_key)
    out = out.merge(host_meta_df, on="Host_id", how="left")
    out = out.sort_values(["PLASMID", "Score"], ascending=[True, False]).copy()

    rows: list[dict[str, object]] = []
    for plasmid, sub in out.groupby("PLASMID", sort=False):
        sub = sub.sort_values("Score", ascending=False).reset_index(drop=True)
        if sub.empty:
            continue
        top1_score = float(sub.loc[0, "Score"])
        cutoff_score = top1_score - float(score_delta)
        selected = sub.loc[sub["Score"] >= cutoff_score].copy()
        if selected.empty:
            selected = sub.head(1).copy()
        selected["Set_rank"] = range(1, len(selected) + 1)
        selected["Predicted_host_top1_score"] = top1_score
        selected["Predicted_host_set_cutoff_score"] = cutoff_score
        selected["Predicted_host_score_delta"] = float(score_delta)
        rows.extend(selected[
            [
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
        ].to_dict("records"))

    if not rows:
        return empty_hostrange_set_df()
    set_df = pd.DataFrame(rows)
    return reorder_columns(set_df, HOSTRANGE_SET_COLUMNS, keep_extra=False)


def build_taxonomy_summary(score_df: pd.DataFrame, model_dir: str, tax_level: str = "species") -> pd.DataFrame:
    if score_df is None or score_df.empty:
        return empty_hostrange_taxonomy_df()
    meta = load_host_metadata(model_dir)
    work = score_df[["PLASMID", "Host_id", "Score"]].copy()
    work["PLASMID"] = work["PLASMID"].map(_norm_plasmid_key)
    return aggregate_scores_by_taxonomy(work, meta, tax_level=tax_level)


def classify_host_range(
    set_df: pd.DataFrame,
    range_threshold_narrow: float | None = None,
    range_threshold_broad: float | None = None,
) -> pd.DataFrame:
    if set_df is None or set_df.empty:
        return pd.DataFrame(columns=["PLASMID", "Predicted_host_range_class"])
    counts = (
        set_df.groupby("PLASMID", as_index=False)
        .size()
        .rename(columns={"size": "Predicted_host_set_size"})
    )
    narrow_cut = 1 if range_threshold_narrow is None else float(range_threshold_narrow)
    broad_cut = 5 if range_threshold_broad is None else float(range_threshold_broad)

    def _label(n: float) -> str:
        if n <= narrow_cut:
            return "narrow"
        if n >= broad_cut:
            return "broad"
        return "intermediate"

    counts["Predicted_host_range_class"] = counts["Predicted_host_set_size"].apply(_label)
    return counts[["PLASMID", "Predicted_host_range_class"]]


def build_hostrange_summary(top_hits_df: pd.DataFrame, set_df: pd.DataFrame, range_df: pd.DataFrame) -> pd.DataFrame:
    if top_hits_df is None or top_hits_df.empty:
        return empty_hostrange_df()

    top1 = (
        top_hits_df.sort_values(["PLASMID", "Rank"])
        .groupby("PLASMID", as_index=False)
        .first()[["PLASMID", "Host_species"]]
        .rename(columns={"Host_species": "Predicted_host_top1"})
    )

    if set_df is None or set_df.empty:
        set_summary = pd.DataFrame(columns=[
            "PLASMID",
            "Predicted_host_set",
            "Predicted_host_set_size",
            "Predicted_host_set_cutoff_score",
            "Predicted_host_score_delta",
            "Predicted_genus_set",
            "Predicted_genus_set_size",
        ])
    else:
        set_work = set_df.copy()
        set_work["PLASMID"] = set_work["PLASMID"].map(_norm_plasmid_key)
        set_summary = (
            set_work.sort_values(["PLASMID", "Set_rank"])
            .groupby("PLASMID", as_index=False)
            .agg(
                Predicted_host_set=("Host_species", lambda x: _first_unique_join(x.tolist())),
                Predicted_host_set_size=("Host_species", lambda x: int(pd.Series(x).astype(str).replace("", pd.NA).dropna().shape[0])),
                Predicted_host_set_cutoff_score=("Predicted_host_set_cutoff_score", "first"),
                Predicted_host_score_delta=("Predicted_host_score_delta", "first"),
                Predicted_genus_set=("Host_genus", lambda x: _first_unique_join(x.tolist())),
                Predicted_genus_set_size=("Host_genus", lambda x: int(pd.Series([str(v).strip() for v in x if str(v).strip()]).drop_duplicates().shape[0])),
            )
        )

    # The model score is reported directly from the top-ranked host.  We do not
    # derive an additional high/medium/low confidence label because the previous
    # heuristic mixed score separation with predicted host-set breadth and was
    # not independently calibrated.
    top_scores = (
        top_hits_df.sort_values(["PLASMID", "Rank"])
        .groupby("PLASMID", as_index=False)
        .first()[["PLASMID", "Score"]]
        .rename(columns={"Score": "Predicted_host_top1_score"})
    )

    summary = top1.merge(top_scores, on="PLASMID", how="left")
    summary = summary.merge(set_summary, on="PLASMID", how="left")
    summary = summary.merge(range_df, on="PLASMID", how="left")
    return reorder_columns(summary, HOSTRANGE_COLUMNS, keep_extra=False)


def postprocess_host_scores(
    score_df: pd.DataFrame,
    model_dir: str,
    topk: int = DEFAULT_TOPK,
    score_delta: float = DEFAULT_SCORE_DELTA,
    tax_level: str = "species",
    range_threshold_narrow: float | None = None,
    range_threshold_broad: float | None = None,
) -> dict[str, pd.DataFrame]:
    if score_df is None or score_df.empty:
        return {
            "top_hits": empty_hostrange_top_hits_df(),
            "set_df": empty_hostrange_set_df(),
            "taxonomy": empty_hostrange_taxonomy_df(),
            "hostrange": empty_hostrange_df(),
        }

    host_meta = load_host_metadata(model_dir)
    top_hits = build_top_hits(score_df, host_meta, topk=topk)
    set_df = build_prediction_sets(score_df, host_meta, score_delta=score_delta)
    taxonomy = build_taxonomy_summary(score_df, model_dir, tax_level=tax_level)
    range_df = classify_host_range(
        set_df,
        range_threshold_narrow=range_threshold_narrow,
        range_threshold_broad=range_threshold_broad,
    )
    hostrange = build_hostrange_summary(top_hits, set_df, range_df)
    return {
        "top_hits": top_hits,
        "set_df": set_df,
        "taxonomy": taxonomy,
        "hostrange": hostrange,
    }
