"""Taxonomy helpers for host-range postprocessing."""

from __future__ import annotations

from pathlib import Path

import joblib
import pandas as pd

from .validate import validate_model_dir

EXPECTED_HOST_META_COLUMNS = ["Host_id", "Host_species", "Host_genus", "Host_family"]


def _fallback_from_vocab(host_vocab_path: str) -> pd.DataFrame:
    obj = joblib.load(host_vocab_path)
    if isinstance(obj, dict) and "id2host_key" in obj:
        id2host = obj["id2host_key"]
    else:
        raise ValueError("Unsupported host vocab format")
    rows = []
    for hid, key in sorted(id2host.items(), key=lambda x: int(x[0]) if str(x[0]).isdigit() else x[0]):
        species = str(key)
        genus = species.split("_")[0] if species else ""
        rows.append({"Host_id": int(hid), "Host_species": species, "Host_genus": genus, "Host_family": ""})
    return pd.DataFrame(rows)


def load_host_metadata(model_dir: str | Path) -> pd.DataFrame:
    files = validate_model_dir(str(model_dir) if model_dir is not None else None)
    path = files.get("host_metadata_path")
    if path and Path(path).is_file():
        df = pd.read_csv(path, sep="\t")
    else:
        df = _fallback_from_vocab(files["host_vocab_path"])
    for col in EXPECTED_HOST_META_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    df["Host_id"] = pd.to_numeric(df["Host_id"], errors="coerce").astype("Int64")
    return normalize_taxonomy_strings(df)


def normalize_taxonomy_strings(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in ("Host_species", "Host_genus", "Host_family"):
        if col in out.columns:
            out[col] = out[col].fillna("").astype(str)
    return out


def aggregate_scores_by_taxonomy(score_df: pd.DataFrame, host_meta_df: pd.DataFrame, tax_level: str) -> pd.DataFrame:
    if tax_level not in {"species", "genus", "family"}:
        raise ValueError(f"Unsupported tax_level: {tax_level}")
    merged = score_df.merge(host_meta_df, on="Host_id", how="left")
    tax_col = {"species": "Host_species", "genus": "Host_genus", "family": "Host_family"}[tax_level]
    agg = (
        merged.groupby(["PLASMID", tax_col], dropna=False)
        .agg(Aggregated_score=("Score", "sum"), N_hosts=("Host_id", "nunique"))
        .reset_index()
        .rename(columns={tax_col: "Taxon"})
    )
    agg["Tax_level"] = tax_level
    agg = agg.sort_values(["PLASMID", "Aggregated_score"], ascending=[True, False])
    agg["Rank"] = agg.groupby("PLASMID").cumcount() + 1
    return agg[["PLASMID", "Tax_level", "Taxon", "Aggregated_score", "Rank", "N_hosts"]]


def get_top_taxon_per_plasmid(taxonomy_df: pd.DataFrame) -> pd.DataFrame:
    if taxonomy_df.empty:
        return taxonomy_df.copy()
    return taxonomy_df.sort_values(["PLASMID", "Rank"]).groupby("PLASMID", as_index=False).first()
