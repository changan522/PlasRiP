"""External reference-database interface for PlasRiP v0.2+.

PlasRiP software and the PTU/PGG reference database are versioned and
released separately.  This module is the single source of truth for the
schema-2.0 database contract and for database integrity checks used by the
CLI, pipeline, and PTU assignment code.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date
import os
from pathlib import Path
import re
import subprocess
from typing import Dict, Iterable, List, Tuple

import pandas as pd

from .config import (
    DATABASE_MANIFEST_FILENAME,
    MASH_INDEX_SUBDIR,
    PGG_METADATA_FILENAME,
    PLASMID_REFERENCE_INDEX_FILENAME,
    PTU_METADATA_FILENAME,
)
from .exceptions import PlasRipConfigError, PlasRipDataError


# ---------------------------------------------------------------------------
# Database schema 2.0
# ---------------------------------------------------------------------------

INDEX_REQUIRED_COLUMNS = [
    "mash_reference_id",
    "plasmid_id",
    "accession",
    "length_bp",
    "PTU_ID",
    "PTU_assignment_status",
    "is_PTU_representative",
]

PTU_REQUIRED_COLUMNS = [
    "PTU_ID",
    "representative_plasmid",
    "n_plasmids",
    "PGG_ID",
    "PGG_assignment_status",
    "mean_length_bp",
    "median_length_bp",
    "min_length_bp",
    "max_length_bp",
    "replicon_annotated_fraction",
    "n_replicon_types",
    "replicon_types",
    "dominant_replicon",
    "dominant_replicon_fraction",
    "n_host_species",
    "observed_host_species",
    "top_host_species",
    "top_host_species_count",
    "top_host_species_coverage",
    "n_host_genera",
    "observed_host_genera",
    "top_host_genus",
    "top_host_genus_count",
    "top_host_genus_coverage",
    "n_host_families",
    "top_host_family",
    "top_host_family_coverage",
    "n_host_orders",
    "top_host_order",
    "top_host_order_coverage",
    "n_host_classes",
    "top_host_class",
    "top_host_class_coverage",
    "n_host_phyla",
    "top_host_phylum",
    "top_host_phylum_count",
    "top_host_phylum_coverage",
    "observed_host_range_level",
    "observed_host_range_rank",
    "dominant_taxon",
    "dominant_taxon_coverage",
    "n_CONJ",
    "CONJ_fraction",
    "n_MOB",
    "MOB_fraction",
    "n_NT",
    "NT_fraction",
    "dominant_mobility",
    "n_ARG_positive_plasmids",
    "ARG_positive_fraction",
    "mean_ARG_count_per_plasmid",
    "n_unique_ARGs",
    "unique_ARGs",
    "n_unique_ARG_classes",
    "unique_ARG_classes",
    "n_virulence_positive_plasmids",
    "virulence_positive_fraction",
    "mean_virulence_gene_count_per_plasmid",
    "n_unique_virulence_genes",
    "unique_virulence_genes",
]

PGG_REQUIRED_COLUMNS = [
    "PGG_ID",
    "n_PTUs",
    "n_plasmids",
    "anchor_PTU",
    "anchor_representative_plasmid",
    "minimum_internal_support",
    "merge_with_parent_support",
    "effective_persistence_width",
    "persistence_volume",
    "anchor_mean_within_group_Mash",
]

MANIFEST_REQUIRED_KEYS = [
    "database_name",
    "database_version",
    "schema_version",
    "release_date",
    "mash_index_directory",
    "total_reference_plasmids",
    "stable_PTU_assigned_plasmids",
    "unassigned_no_stable_PTU_plasmids",
    "PTU_count",
    "PGG_assigned_PTUs",
    "PGG_count",
    "PTU_ID_format",
    "PGG_ID_format",
    "mash_version",
    "mash_hash_function",
    "mash_seed",
    "mash_k",
    "mash_alphabet",
    "mash_target_min_hashes_per_sketch",
    "mash_sketches_per_file",
    "mash_reference_id_pattern",
    "query_PTU_Mash_cutoff",
]

SUPPORTED_SCHEMA_VERSION = "2.0"
EXPECTED_PTU_ID_EXAMPLE = "PTU000001"
EXPECTED_PGG_ID_EXAMPLE = "PGG-000001"
EXPECTED_MASH_HASH_FUNCTION = "MurmurHash3_x64_128"
EXPECTED_MASH_ALPHABET = "ACGT (canonical)"
EXPECTED_SKETCHES_PER_FILE = 1

PTU_ID_RE = re.compile(r"^PTU\d{6}$")
PGG_ID_RE = re.compile(r"^PGG-\d{6}$")

INDEX_ASSIGNED_STATUS = "ASSIGNED_STABLE_PTU"
INDEX_UNASSIGNED_STATUS = "UNASSIGNED_NO_STABLE_PTU"
PGG_ASSIGNED_STATUS = "ASSIGNED_STABLE_PGG"
PGG_UNASSIGNED_STATUSES = {
    "UNASSIGNED_OUTSIDE_CONSENSUS_SPACE",
    "UNASSIGNED_NOT_ELIGIBLE",
    "UNASSIGNED_NO_STABLE_GROUP",
}

CANONICAL_BOOL_VALUES = {"TRUE", "FALSE"}

PTU_FRACTION_COLUMNS = [
    "replicon_annotated_fraction",
    "dominant_replicon_fraction",
    "top_host_species_coverage",
    "top_host_genus_coverage",
    "top_host_family_coverage",
    "top_host_order_coverage",
    "top_host_class_coverage",
    "top_host_phylum_coverage",
    "dominant_taxon_coverage",
    "CONJ_fraction",
    "MOB_fraction",
    "NT_fraction",
    "ARG_positive_fraction",
    "virulence_positive_fraction",
]

HOST_RANGE_LEVEL_TO_RANK = {
    "1": "Species",
    "2": "Genus",
    "3": "Family",
    "4": "Order",
    "5": "Class",
    "6": "Phylum",
    "7": "Cross-phylum",
}


@dataclass(frozen=True)
class MashSketchInfo:
    """Mash sketch header information parsed from ``mash info``."""

    hash_function: str
    seed: int
    k: int
    alphabet: str
    sketch_size: int
    sketches: int


@dataclass(frozen=True)
class PlasRipReferenceDatabase:
    """Loaded PlasRiP schema-2.0 reference database."""

    root: Path
    mash_index_dir: Path
    plasmid_index: pd.DataFrame
    ptu_metadata: pd.DataFrame
    pgg_metadata: pd.DataFrame
    manifest: Dict[str, str]

    @property
    def database_version(self) -> str:
        return self.manifest["database_version"]

    @property
    def schema_version(self) -> str:
        return self.manifest["schema_version"]

    @property
    def mash_version(self) -> str:
        return self.manifest["mash_version"]

    @property
    def mash_hash_function(self) -> str:
        return self.manifest["mash_hash_function"]

    @property
    def mash_seed(self) -> int:
        return int(self.manifest["mash_seed"])

    @property
    def mash_k(self) -> int:
        return int(self.manifest["mash_k"])

    @property
    def mash_sketch_size(self) -> int:
        return int(self.manifest["mash_target_min_hashes_per_sketch"])

    @property
    def mash_alphabet(self) -> str:
        return self.manifest["mash_alphabet"]

    @property
    def mash_reference_id_pattern(self) -> str:
        return self.manifest["mash_reference_id_pattern"]

    @property
    def ptu_threshold(self) -> float:
        return float(self.manifest["query_PTU_Mash_cutoff"])


# ---------------------------------------------------------------------------
# File loading and primitive validation helpers
# ---------------------------------------------------------------------------


def resolve_database_dir(database_dir: str | os.PathLike[str] | None = None) -> Path:
    """Resolve the external reference-database directory.

    Resolution priority:
      1. explicit ``database_dir``
      2. ``PLASRIP_DATABASE`` environment variable

    No PTU/PGG reference database is bundled inside the Python package.
    """
    if database_dir:
        return Path(database_dir).expanduser().resolve()

    env_value = os.environ.get("PLASRIP_DATABASE")
    if env_value:
        return Path(env_value).expanduser().resolve()

    raise PlasRipConfigError(
        "PlasRiP reference database was not specified. Use --database <DIR> "
        "or set PLASRIP_DATABASE=<DIR>."
    )


def _read_tsv(path: Path, required_columns: Iterable[str], label: str) -> pd.DataFrame:
    if not path.is_file():
        raise PlasRipDataError(f"Missing {label}: {path}")
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            encoding="utf-8-sig",
            dtype=str,
            keep_default_na=False,
        )
    except Exception as exc:  # pragma: no cover - pandas wording depends on version
        raise PlasRipDataError(f"Could not read {label}: {path}: {exc}") from exc

    df.columns = df.columns.str.replace("\ufeff", "", regex=False).str.strip()
    missing = [c for c in required_columns if c not in df.columns]
    if missing:
        raise PlasRipDataError(
            f"{label} is missing required schema-2.0 columns {missing}. "
            f"Found: {list(df.columns)}"
        )
    return df


def load_manifest(path: Path) -> Dict[str, str]:
    df = _read_tsv(path, ["key", "value"], "database manifest")
    df["key"] = df["key"].astype(str).str.strip()
    df["value"] = df["value"].astype(str).str.strip()

    if (df["key"] == "").any():
        raise PlasRipDataError("Database manifest contains an empty key")
    if df["key"].duplicated().any():
        dup = df.loc[df["key"].duplicated(), "key"].tolist()
        raise PlasRipDataError(f"Duplicate keys in database manifest: {dup}")

    manifest = dict(zip(df["key"], df["value"]))
    missing = [k for k in MANIFEST_REQUIRED_KEYS if not manifest.get(k)]
    if missing:
        raise PlasRipDataError(f"Database manifest is missing required keys: {missing}")
    return manifest


def _parse_nonnegative_int(value: str, label: str) -> int:
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise PlasRipDataError(f"{label} must be an integer, got {value!r}") from exc
    if parsed < 0:
        raise PlasRipDataError(f"{label} must be non-negative, got {parsed}")
    return parsed


def _parse_positive_int(value: str, label: str) -> int:
    parsed = _parse_nonnegative_int(value, label)
    if parsed <= 0:
        raise PlasRipDataError(f"{label} must be > 0, got {parsed}")
    return parsed


def _parse_float(value: str, label: str) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError) as exc:
        raise PlasRipDataError(f"{label} must be numeric, got {value!r}") from exc
    return parsed


def _parse_fraction(value: str, label: str, *, allow_blank: bool = False) -> float | None:
    if allow_blank and str(value) == "":
        return None
    parsed = _parse_float(value, label)
    if not 0.0 <= parsed <= 1.0:
        raise PlasRipDataError(f"{label} must be within [0, 1], got {parsed}")
    return parsed


def _validate_id_format(values: Iterable[str], pattern: re.Pattern[str], label: str) -> None:
    bad = [str(v) for v in values if str(v) and pattern.fullmatch(str(v)) is None]
    if bad:
        raise PlasRipDataError(f"Invalid {label} format (showing up to 10): {bad[:10]}")


def _validate_nonempty(df: pd.DataFrame, columns: Iterable[str], label: str) -> None:
    for column in columns:
        bad_mask = df[column].astype(str).str.strip() == ""
        if bad_mask.any():
            examples = df.loc[bad_mask].head(10).to_dict(orient="records")
            raise PlasRipDataError(
                f"{label} contains empty values in required field {column!r}. "
                f"Examples: {examples}"
            )


def _validate_canonical_bool_column(df: pd.DataFrame, column: str, label: str) -> None:
    values = df[column].astype(str).str.strip().str.upper()
    bad = ~values.isin(CANONICAL_BOOL_VALUES)
    if bad.any():
        examples = df.loc[bad, [column]].head(10).to_dict(orient="records")
        raise PlasRipDataError(
            f"{label} field {column!r} must use canonical TRUE/FALSE values. "
            f"Examples: {examples}"
        )


def _as_bool_series(s: pd.Series) -> pd.Series:
    """Convert a previously validated canonical TRUE/FALSE field to bool."""
    return s.astype(str).str.strip().str.upper() == "TRUE"


def _resolve_mash_index_dir(root: Path, manifest: Dict[str, str]) -> Path:
    raw = manifest["mash_index_directory"].strip()
    rel = Path(raw)
    if rel.is_absolute():
        raise PlasRipDataError(
            "mash_index_directory must be relative to the database root, "
            f"got absolute path: {raw!r}"
        )

    root_resolved = root.resolve()
    mash_dir = (root_resolved / rel).resolve()
    try:
        mash_dir.relative_to(root_resolved)
    except ValueError as exc:
        raise PlasRipDataError(
            "mash_index_directory escapes the database root: "
            f"{raw!r} -> {mash_dir}"
        ) from exc
    return mash_dir


def _validate_manifest_semantics(manifest: Dict[str, str]) -> None:
    if manifest["schema_version"] != SUPPORTED_SCHEMA_VERSION:
        raise PlasRipDataError(
            f"Unsupported database schema_version={manifest['schema_version']!r}; "
            f"PlasRiP v0.2 supports schema {SUPPORTED_SCHEMA_VERSION}."
        )

    try:
        date.fromisoformat(manifest["release_date"])
    except ValueError as exc:
        raise PlasRipDataError(
            f"Manifest release_date must use ISO YYYY-MM-DD format, got {manifest['release_date']!r}"
        ) from exc

    if manifest["PTU_ID_format"] != EXPECTED_PTU_ID_EXAMPLE:
        raise PlasRipDataError(
            f"Unsupported PTU_ID_format={manifest['PTU_ID_format']!r}; "
            f"schema 2.0 requires {EXPECTED_PTU_ID_EXAMPLE!r}."
        )
    if manifest["PGG_ID_format"] != EXPECTED_PGG_ID_EXAMPLE:
        raise PlasRipDataError(
            f"Unsupported PGG_ID_format={manifest['PGG_ID_format']!r}; "
            f"schema 2.0 requires {EXPECTED_PGG_ID_EXAMPLE!r}."
        )

    if manifest["mash_hash_function"] != EXPECTED_MASH_HASH_FUNCTION:
        raise PlasRipDataError(
            f"Unsupported Mash hash function {manifest['mash_hash_function']!r}; "
            f"schema 2.0 expects {EXPECTED_MASH_HASH_FUNCTION!r}."
        )
    if manifest["mash_alphabet"] != EXPECTED_MASH_ALPHABET:
        raise PlasRipDataError(
            f"Unsupported Mash alphabet {manifest['mash_alphabet']!r}; "
            f"schema 2.0 expects {EXPECTED_MASH_ALPHABET!r}."
        )

    _parse_nonnegative_int(manifest["mash_seed"], "manifest mash_seed")
    _parse_positive_int(manifest["mash_k"], "manifest mash_k")
    _parse_positive_int(
        manifest["mash_target_min_hashes_per_sketch"],
        "manifest mash_target_min_hashes_per_sketch",
    )
    sketches_per_file = _parse_positive_int(
        manifest["mash_sketches_per_file"], "manifest mash_sketches_per_file"
    )
    if sketches_per_file != EXPECTED_SKETCHES_PER_FILE:
        raise PlasRipDataError(
            f"schema 2.0 requires exactly {EXPECTED_SKETCHES_PER_FILE} sketch per .msh file; "
            f"manifest reports {sketches_per_file}."
        )

    cutoff = _parse_float(manifest["query_PTU_Mash_cutoff"], "manifest query_PTU_Mash_cutoff")
    if not 0.0 < cutoff < 1.0:
        raise PlasRipDataError(
            f"manifest query_PTU_Mash_cutoff must be in (0, 1), got {cutoff}"
        )

    if not manifest["mash_reference_id_pattern"].strip():
        raise PlasRipDataError("manifest mash_reference_id_pattern must not be empty")


# ---------------------------------------------------------------------------
# Cross-table database validation
# ---------------------------------------------------------------------------


def validate_database_dir(
    database_dir: str | os.PathLike[str] | None,
    *,
    check_sketches: bool = False,
    check_all_sketches: bool = False,
) -> Dict[str, object]:
    """Validate a PlasRiP schema-2.0 reference database.

    Default validation checks file presence, schema, manifest semantics, and
    cross-table integrity without walking the Mash index.

    ``check_sketches=True`` additionally verifies every PTU representative
    sketch. ``check_all_sketches=True`` verifies all reference-plasmid sketches
    and is intended for full release QC.
    """
    root = resolve_database_dir(database_dir)
    if not root.is_dir():
        raise PlasRipConfigError(f"Reference database directory not found: {root}")

    manifest_path = root / DATABASE_MANIFEST_FILENAME
    index_path = root / PLASMID_REFERENCE_INDEX_FILENAME
    ptu_path = root / PTU_METADATA_FILENAME
    pgg_path = root / PGG_METADATA_FILENAME

    manifest = load_manifest(manifest_path)
    _validate_manifest_semantics(manifest)

    mash_dir = _resolve_mash_index_dir(root, manifest)
    if not mash_dir.is_dir():
        raise PlasRipDataError(f"Mash index directory not found: {mash_dir}")

    index = _read_tsv(index_path, INDEX_REQUIRED_COLUMNS, "plasmid reference index")
    ptu = _read_tsv(ptu_path, PTU_REQUIRED_COLUMNS, "PTU metadata")
    pgg = _read_tsv(pgg_path, PGG_REQUIRED_COLUMNS, "PGG metadata")

    # ----- nonempty key fields, uniqueness, and identifier format -----
    _validate_nonempty(
        index,
        ["mash_reference_id", "plasmid_id", "accession", "length_bp", "PTU_assignment_status", "is_PTU_representative"],
        "plasmid reference index",
    )
    _validate_nonempty(
        ptu,
        ["PTU_ID", "representative_plasmid", "n_plasmids", "PGG_assignment_status"],
        "PTU metadata",
    )
    _validate_nonempty(
        pgg,
        ["PGG_ID", "n_PTUs", "n_plasmids", "anchor_PTU", "anchor_representative_plasmid"],
        "PGG metadata",
    )

    if index["plasmid_id"].duplicated().any():
        raise PlasRipDataError("plasmid_reference_index.tsv contains duplicate plasmid_id values")
    if index["mash_reference_id"].duplicated().any():
        raise PlasRipDataError("plasmid_reference_index.tsv contains duplicate mash_reference_id values")
    if ptu["PTU_ID"].duplicated().any():
        raise PlasRipDataError("PTU_metadata.tsv contains duplicate PTU_ID values")
    if pgg["PGG_ID"].duplicated().any():
        raise PlasRipDataError("PGG_metadata.tsv contains duplicate PGG_ID values")

    _validate_id_format(ptu["PTU_ID"], PTU_ID_RE, "PTU_ID")
    _validate_id_format(pgg["PGG_ID"], PGG_ID_RE, "PGG_ID")
    _validate_canonical_bool_column(index, "is_PTU_representative", "plasmid reference index")

    # Positive plasmid lengths are part of the runtime reference-index contract.
    invalid_lengths = []
    for row in index[["plasmid_id", "length_bp"]].itertuples(index=False):
        try:
            _parse_positive_int(str(row.length_bp), f"plasmid {row.plasmid_id} length_bp")
        except PlasRipDataError:
            invalid_lengths.append((str(row.plasmid_id), str(row.length_bp)))
    if invalid_lengths:
        raise PlasRipDataError(
            "plasmid_reference_index.tsv contains invalid length_bp values "
            f"(showing up to 10): {invalid_lengths[:10]}"
        )

    # ----- manifest counts -----
    expected_reference = _parse_nonnegative_int(
        manifest["total_reference_plasmids"], "manifest total_reference_plasmids"
    )
    expected_stable = _parse_nonnegative_int(
        manifest["stable_PTU_assigned_plasmids"], "manifest stable_PTU_assigned_plasmids"
    )
    expected_unassigned = _parse_nonnegative_int(
        manifest["unassigned_no_stable_PTU_plasmids"],
        "manifest unassigned_no_stable_PTU_plasmids",
    )
    expected_ptu = _parse_nonnegative_int(manifest["PTU_count"], "manifest PTU_count")
    expected_pgg_assigned = _parse_nonnegative_int(
        manifest["PGG_assigned_PTUs"], "manifest PGG_assigned_PTUs"
    )
    expected_pgg = _parse_nonnegative_int(manifest["PGG_count"], "manifest PGG_count")

    if len(index) != expected_reference:
        raise PlasRipDataError(
            f"Reference-plasmid count mismatch: manifest={expected_reference}, "
            f"plasmid_reference_index.tsv={len(index)}"
        )
    if len(ptu) != expected_ptu:
        raise PlasRipDataError(
            f"PTU count mismatch: manifest={expected_ptu}, PTU_metadata.tsv={len(ptu)}"
        )
    if len(pgg) != expected_pgg:
        raise PlasRipDataError(
            f"PGG count mismatch: manifest={expected_pgg}, PGG_metadata.tsv={len(pgg)}"
        )
    if expected_stable + expected_unassigned != expected_reference:
        raise PlasRipDataError(
            "Manifest plasmid counts are inconsistent: "
            f"stable ({expected_stable}) + unassigned ({expected_unassigned}) "
            f"!= total ({expected_reference})."
        )

    # ----- plasmid -> PTU integrity -----
    index_ptu = index["PTU_ID"].astype(str).str.strip()
    index_status = index["PTU_assignment_status"].astype(str).str.strip()
    assigned_mask = index_ptu != ""
    stable_mask = index_status == INDEX_ASSIGNED_STATUS
    unassigned_mask = index_status == INDEX_UNASSIGNED_STATUS

    invalid_index_status = index.loc[
        ~index_status.isin({INDEX_ASSIGNED_STATUS, INDEX_UNASSIGNED_STATUS})
    ]
    if not invalid_index_status.empty:
        bad = invalid_index_status[["plasmid_id", "PTU_assignment_status"]].head(10).to_dict(orient="records")
        raise PlasRipDataError(f"Invalid PTU_assignment_status values: {bad}")
    if not (assigned_mask == stable_mask).all():
        bad = index.loc[
            assigned_mask != stable_mask,
            ["plasmid_id", "PTU_ID", "PTU_assignment_status"],
        ].head(10)
        raise PlasRipDataError(
            "PTU_ID/PTU_assignment_status inconsistency in plasmid reference index. "
            f"Examples: {bad.to_dict(orient='records')}"
        )
    if not ((~assigned_mask) == unassigned_mask).all():
        bad = index.loc[
            (~assigned_mask) != unassigned_mask,
            ["plasmid_id", "PTU_ID", "PTU_assignment_status"],
        ].head(10)
        raise PlasRipDataError(
            "Unassigned plasmid status is inconsistent with empty PTU_ID. "
            f"Examples: {bad.to_dict(orient='records')}"
        )
    if int(stable_mask.sum()) != expected_stable:
        raise PlasRipDataError(
            f"Stable PTU-assigned plasmid count mismatch: manifest={expected_stable}, "
            f"index={int(stable_mask.sum())}"
        )
    if int(unassigned_mask.sum()) != expected_unassigned:
        raise PlasRipDataError(
            f"Unassigned plasmid count mismatch: manifest={expected_unassigned}, "
            f"index={int(unassigned_mask.sum())}"
        )

    ptu_ids = set(ptu["PTU_ID"].astype(str))
    _validate_id_format(index.loc[assigned_mask, "PTU_ID"], PTU_ID_RE, "index PTU_ID")
    bad_ptu = sorted(set(index.loc[assigned_mask, "PTU_ID"].astype(str)) - ptu_ids)
    if bad_ptu:
        raise PlasRipDataError(
            "Reference index contains PTU IDs absent from PTU metadata "
            f"(showing up to 10): {bad_ptu[:10]}"
        )

    # Every PTU must have the number of stable members reported in metadata.
    member_counts = index.loc[stable_mask].groupby("PTU_ID").size().to_dict()
    metadata_member_counts: Dict[str, int] = {}
    for row in ptu[["PTU_ID", "n_plasmids"]].itertuples(index=False):
        metadata_member_counts[str(row.PTU_ID)] = _parse_positive_int(
            str(row.n_plasmids), f"PTU {row.PTU_ID} n_plasmids"
        )

    missing_ptu_members = sorted(ptu_ids - set(member_counts))
    if missing_ptu_members:
        raise PlasRipDataError(
            f"PTUs without stable member plasmids (showing up to 10): {missing_ptu_members[:10]}"
        )
    member_mismatch = [
        (ptu_id, metadata_member_counts[ptu_id], int(member_counts.get(ptu_id, 0)))
        for ptu_id in sorted(ptu_ids)
        if metadata_member_counts[ptu_id] != int(member_counts.get(ptu_id, 0))
    ]
    if member_mismatch:
        raise PlasRipDataError(
            "PTU member-count mismatch between PTU_metadata.tsv and plasmid_reference_index.tsv "
            f"(showing up to 10): {member_mismatch[:10]}"
        )

    # ----- PTU metadata numeric sanity -----
    for row in ptu[["PTU_ID", "min_length_bp", "median_length_bp", "mean_length_bp", "max_length_bp"]].itertuples(index=False):
        ptu_id = str(row.PTU_ID)
        min_len = _parse_float(str(row.min_length_bp), f"PTU {ptu_id} min_length_bp")
        median_len = _parse_float(str(row.median_length_bp), f"PTU {ptu_id} median_length_bp")
        mean_len = _parse_float(str(row.mean_length_bp), f"PTU {ptu_id} mean_length_bp")
        max_len = _parse_float(str(row.max_length_bp), f"PTU {ptu_id} max_length_bp")
        if min_len <= 0 or not (min_len <= median_len <= max_len) or not (min_len <= mean_len <= max_len):
            raise PlasRipDataError(
                f"Invalid PTU length summary for {ptu_id}: "
                f"min={min_len}, median={median_len}, mean={mean_len}, max={max_len}"
            )

    optional_fraction_count_columns = {
        "dominant_replicon_fraction": "n_replicon_types",
        "top_host_species_coverage": "n_host_species",
        "top_host_genus_coverage": "n_host_genera",
        "top_host_family_coverage": "n_host_families",
        "top_host_order_coverage": "n_host_orders",
        "top_host_class_coverage": "n_host_classes",
        "top_host_phylum_coverage": "n_host_phyla",
    }
    for column in PTU_FRACTION_COLUMNS:
        bad = []
        count_column = optional_fraction_count_columns.get(column)
        columns = ["PTU_ID", column] + ([count_column] if count_column else [])
        for row in ptu[columns].itertuples(index=False):
            value = getattr(row, column)
            if value == "" and count_column:
                n_available = _parse_nonnegative_int(
                    str(getattr(row, count_column)),
                    f"PTU {row.PTU_ID} {count_column}",
                )
                if n_available == 0:
                    continue
            try:
                _parse_fraction(str(value), f"PTU {row.PTU_ID} {column}")
            except PlasRipDataError:
                bad.append((str(row.PTU_ID), str(value)))
        if bad:
            raise PlasRipDataError(
                f"PTU_metadata.tsv contains invalid {column} values (showing up to 10): {bad[:10]}"
            )

    # Mobility is exhaustive in the schema-2.0 source metadata.
    bad_mobility = []
    for row in ptu[["PTU_ID", "n_plasmids", "n_CONJ", "n_MOB", "n_NT"]].itertuples(index=False):
        n = metadata_member_counts[str(row.PTU_ID)]
        counts = [
            _parse_nonnegative_int(str(row.n_CONJ), f"PTU {row.PTU_ID} n_CONJ"),
            _parse_nonnegative_int(str(row.n_MOB), f"PTU {row.PTU_ID} n_MOB"),
            _parse_nonnegative_int(str(row.n_NT), f"PTU {row.PTU_ID} n_NT"),
        ]
        if sum(counts) != n:
            bad_mobility.append((str(row.PTU_ID), n, *counts))
    if bad_mobility:
        raise PlasRipDataError(
            "PTU mobility counts do not sum to n_plasmids (showing up to 10): "
            f"{bad_mobility[:10]}"
        )

    # Host-range rank/level fields must agree.
    bad_host_rank = []
    for row in ptu[["PTU_ID", "observed_host_range_level", "observed_host_range_rank"]].itertuples(index=False):
        level = str(row.observed_host_range_level)
        rank = str(row.observed_host_range_rank)
        if HOST_RANGE_LEVEL_TO_RANK.get(level) != rank:
            bad_host_rank.append((str(row.PTU_ID), level, rank))
    if bad_host_rank:
        raise PlasRipDataError(
            "PTU observed_host_range_level/rank mismatch (showing up to 10): "
            f"{bad_host_rank[:10]}"
        )

    # ----- exactly one representative per PTU, and metadata must point to it -----
    rep_mask = _as_bool_series(index["is_PTU_representative"])
    representative = index.loc[rep_mask].copy()
    if len(representative) != len(ptu):
        raise PlasRipDataError(
            f"Representative count mismatch: representatives={len(representative)}, PTUs={len(ptu)}"
        )
    if (representative["PTU_ID"] == "").any():
        bad = representative.loc[representative["PTU_ID"] == "", "plasmid_id"].head(10).tolist()
        raise PlasRipDataError(f"PTU representatives without a PTU assignment: {bad}")
    rep_counts = representative.groupby("PTU_ID").size()
    bad_rep_counts = rep_counts[rep_counts != 1]
    if len(rep_counts) != len(ptu) or not bad_rep_counts.empty:
        raise PlasRipDataError(
            "Each PTU must have exactly one representative plasmid. "
            f"Non-unit representative counts: {bad_rep_counts.head(10).to_dict()}"
        )
    rep_by_ptu = representative.set_index("PTU_ID")["plasmid_id"].to_dict()
    bad_representative = []
    for row in ptu[["PTU_ID", "representative_plasmid"]].itertuples(index=False):
        if rep_by_ptu.get(str(row.PTU_ID)) != str(row.representative_plasmid):
            bad_representative.append(
                (str(row.PTU_ID), str(row.representative_plasmid), rep_by_ptu.get(str(row.PTU_ID)))
            )
    if bad_representative:
        raise PlasRipDataError(
            "PTU representative mismatch between PTU_metadata.tsv and plasmid_reference_index.tsv "
            f"(showing up to 10): {bad_representative[:10]}"
        )

    # ----- PTU -> PGG integrity -----
    ptu_pgg_id = ptu["PGG_ID"].astype(str).str.strip()
    ptu_pgg_status = ptu["PGG_assignment_status"].astype(str).str.strip()
    valid_pgg_statuses = {PGG_ASSIGNED_STATUS} | PGG_UNASSIGNED_STATUSES
    invalid_pgg_status = ptu.loc[~ptu_pgg_status.isin(valid_pgg_statuses)]
    if not invalid_pgg_status.empty:
        bad = invalid_pgg_status[["PTU_ID", "PGG_assignment_status"]].head(10).to_dict(orient="records")
        raise PlasRipDataError(f"Invalid PGG_assignment_status values: {bad}")

    pgg_assigned_mask = ptu_pgg_status == PGG_ASSIGNED_STATUS
    if not ((ptu_pgg_id != "") == pgg_assigned_mask).all():
        bad = ptu.loc[
            (ptu_pgg_id != "") != pgg_assigned_mask,
            ["PTU_ID", "PGG_ID", "PGG_assignment_status"],
        ].head(10)
        raise PlasRipDataError(
            "PGG_ID/PGG_assignment_status inconsistency in PTU metadata. "
            f"Examples: {bad.to_dict(orient='records')}"
        )
    if int(pgg_assigned_mask.sum()) != expected_pgg_assigned:
        raise PlasRipDataError(
            f"PGG-assigned PTU count mismatch: manifest={expected_pgg_assigned}, "
            f"PTU metadata={int(pgg_assigned_mask.sum())}"
        )

    pgg_ids = set(pgg["PGG_ID"].astype(str))
    _validate_id_format(ptu.loc[pgg_assigned_mask, "PGG_ID"], PGG_ID_RE, "PTU metadata PGG_ID")
    bad_pgg = sorted(set(ptu.loc[pgg_assigned_mask, "PGG_ID"].astype(str)) - pgg_ids)
    if bad_pgg:
        raise PlasRipDataError(
            f"PTU metadata contains PGG IDs absent from PGG metadata (showing up to 10): {bad_pgg[:10]}"
        )

    ptu_count_by_pgg = ptu.loc[pgg_assigned_mask].groupby("PGG_ID").size().to_dict()
    plasmid_count_by_pgg = (
        ptu.loc[pgg_assigned_mask, ["PTU_ID", "PGG_ID"]]
        .assign(n_plasmids=lambda d: d["PTU_ID"].map(metadata_member_counts))
        .groupby("PGG_ID")["n_plasmids"]
        .sum()
        .to_dict()
    )
    ptu_rep_lookup = ptu.set_index("PTU_ID")["representative_plasmid"].to_dict()

    pgg_mismatch = []
    anchor_mismatch = []
    for row in pgg[
        ["PGG_ID", "n_PTUs", "n_plasmids", "anchor_PTU", "anchor_representative_plasmid"]
    ].itertuples(index=False):
        pgg_id = str(row.PGG_ID)
        n_ptus = _parse_positive_int(str(row.n_PTUs), f"PGG {pgg_id} n_PTUs")
        n_plasmids = _parse_positive_int(str(row.n_plasmids), f"PGG {pgg_id} n_plasmids")
        actual_ptus = int(ptu_count_by_pgg.get(pgg_id, 0))
        actual_plasmids = int(plasmid_count_by_pgg.get(pgg_id, 0))
        if n_ptus != actual_ptus or n_plasmids != actual_plasmids:
            pgg_mismatch.append((pgg_id, n_ptus, actual_ptus, n_plasmids, actual_plasmids))

        anchor_ptu = str(row.anchor_PTU)
        anchor_rep = str(row.anchor_representative_plasmid)
        anchor_ok = (
            anchor_ptu in ptu_ids
            and ptu.loc[ptu["PTU_ID"] == anchor_ptu, "PGG_ID"].iloc[0] == pgg_id
            and ptu_rep_lookup.get(anchor_ptu) == anchor_rep
        )
        if not anchor_ok:
            anchor_mismatch.append((pgg_id, anchor_ptu, anchor_rep))

    if pgg_mismatch:
        raise PlasRipDataError(
            "PGG size mismatch between PGG_metadata.tsv and PTU_metadata.tsv "
            "(PGG, expected_PTUs, actual_PTUs, expected_plasmids, actual_plasmids; up to 10): "
            f"{pgg_mismatch[:10]}"
        )
    if anchor_mismatch:
        raise PlasRipDataError(
            "Invalid PGG anchor PTU/representative mapping (showing up to 10): "
            f"{anchor_mismatch[:10]}"
        )

    # PGG stability metrics are required parts of schema 2.0.
    for row in pgg.itertuples(index=False):
        pgg_id = str(row.PGG_ID)
        support = _parse_fraction(
            str(row.minimum_internal_support), f"PGG {pgg_id} minimum_internal_support"
        )
        parent_support = _parse_fraction(
            str(row.merge_with_parent_support), f"PGG {pgg_id} merge_with_parent_support"
        )
        width = _parse_float(
            str(row.effective_persistence_width), f"PGG {pgg_id} effective_persistence_width"
        )
        volume = _parse_float(str(row.persistence_volume), f"PGG {pgg_id} persistence_volume")
        mash = _parse_fraction(
            str(row.anchor_mean_within_group_Mash), f"PGG {pgg_id} anchor_mean_within_group_Mash"
        )
        if width < 0 or volume < 0:
            raise PlasRipDataError(
                f"PGG {pgg_id} persistence metrics must be non-negative: width={width}, volume={volume}"
            )
        # Variables are assigned to make explicit that parsing is intentional.
        _ = (support, parent_support, mash)

    # ----- optional Mash-index filesystem validation -----
    missing_reps: List[str] = []
    missing_all: List[str] = []
    if check_sketches or check_all_sketches:
        for plasmid_id in representative["plasmid_id"]:
            path = mash_dir / f"{plasmid_id}.msh"
            if not path.is_file():
                missing_reps.append(str(path))
        if missing_reps:
            raise PlasRipDataError(
                f"Missing representative Mash sketches: {len(missing_reps)}. Examples: {missing_reps[:5]}"
            )

    if check_all_sketches:
        for plasmid_id in index["plasmid_id"]:
            path = mash_dir / f"{plasmid_id}.msh"
            if not path.is_file():
                missing_all.append(str(path))
        if missing_all:
            raise PlasRipDataError(
                f"Missing reference Mash sketches: {len(missing_all)}. Examples: {missing_all[:5]}"
            )

    return {
        "root": str(root),
        "database_version": manifest["database_version"],
        "schema_version": manifest["schema_version"],
        "mash_index_dir": str(mash_dir),
        "reference_plasmids": len(index),
        "stable_PTU_assigned_plasmids": int(stable_mask.sum()),
        "unassigned_plasmids": int(unassigned_mask.sum()),
        "PTUs": len(ptu),
        "PGG_assigned_PTUs": int(pgg_assigned_mask.sum()),
        "PGGs": len(pgg),
        "representatives": len(representative),
        "mash_version": manifest["mash_version"],
        "mash_hash_function": manifest["mash_hash_function"],
        "mash_seed": int(manifest["mash_seed"]),
        "mash_k": int(manifest["mash_k"]),
        "mash_sketch_size": int(manifest["mash_target_min_hashes_per_sketch"]),
        "mash_alphabet": manifest["mash_alphabet"],
        "ptu_threshold": float(manifest["query_PTU_Mash_cutoff"]),
        "representative_sketches_checked": bool(check_sketches or check_all_sketches),
        "all_sketches_checked": bool(check_all_sketches),
    }


# ---------------------------------------------------------------------------
# Database loading and runtime lookup helpers
# ---------------------------------------------------------------------------


def load_database(
    database_dir: str | os.PathLike[str] | None,
    *,
    check_sketches: bool = True,
) -> PlasRipReferenceDatabase:
    """Validate and load the external PlasRiP reference database."""
    report = validate_database_dir(database_dir, check_sketches=check_sketches)
    root = Path(str(report["root"]))
    manifest = load_manifest(root / DATABASE_MANIFEST_FILENAME)
    mash_dir = _resolve_mash_index_dir(root, manifest)

    index = _read_tsv(
        root / PLASMID_REFERENCE_INDEX_FILENAME,
        INDEX_REQUIRED_COLUMNS,
        "plasmid reference index",
    )
    ptu = _read_tsv(root / PTU_METADATA_FILENAME, PTU_REQUIRED_COLUMNS, "PTU metadata")
    pgg = _read_tsv(root / PGG_METADATA_FILENAME, PGG_REQUIRED_COLUMNS, "PGG metadata")

    return PlasRipReferenceDatabase(
        root=root,
        mash_index_dir=mash_dir,
        plasmid_index=index,
        ptu_metadata=ptu,
        pgg_metadata=pgg,
        manifest=manifest,
    )


def validate_query_mash_parameters(
    db: PlasRipReferenceDatabase,
    mash_k: int | None,
    mash_s: int | None,
) -> Tuple[int, int]:
    """Return database-compatible query sketch parameters.

    User overrides are accepted only when they exactly match the database.
    This prevents silently computing invalid distances against incompatible
    reference sketches.
    """
    db_k = db.mash_k
    db_s = db.mash_sketch_size

    if mash_k is not None and int(mash_k) != db_k:
        raise PlasRipConfigError(
            f"Mash k-mer size mismatch: database requires k={db_k}, requested k={mash_k}."
        )
    if mash_s is not None and int(mash_s) != db_s:
        raise PlasRipConfigError(
            f"Mash sketch-size mismatch: database requires s={db_s}, requested s={mash_s}."
        )
    return db_k, db_s



_VERSION_RE = re.compile(r"(\d+(?:\.\d+)+)")


def _version_tuple(value: str) -> Tuple[int, ...]:
    match = _VERSION_RE.search(str(value))
    if not match:
        raise PlasRipConfigError(f"Could not parse Mash version from: {value!r}")
    return tuple(int(x) for x in match.group(1).split("."))


def get_mash_runtime_version(mash_cmd: str = "mash") -> str:
    """Return the installed Mash version string.

    The executable is checked explicitly before any PTU assignment starts so
    users receive a clear configuration error instead of a later subprocess
    failure.
    """
    try:
        proc = subprocess.run(
            [mash_cmd, "--version"],
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as exc:
        raise PlasRipConfigError(
            f"Mash executable not found: {mash_cmd!r}. Install Mash or provide the correct executable."
        ) from exc
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or exc.stdout or "").strip()
        raise PlasRipConfigError(
            f"Failed to execute '{mash_cmd} --version'. {detail}".strip()
        ) from exc

    output = (proc.stdout or proc.stderr or "").strip()
    match = _VERSION_RE.search(output)
    if not match:
        raise PlasRipConfigError(
            f"Could not parse Mash version from '{mash_cmd} --version' output: {output!r}"
        )
    return match.group(1)


def validate_mash_runtime(
    db: PlasRipReferenceDatabase,
    mash_cmd: str = "mash",
) -> str:
    """Validate the runtime Mash release against the database manifest.

    The database records a release at major/minor precision (currently 2.3).
    A runtime patch release such as 2.3.1 is therefore accepted, while a
    different major/minor release is rejected for reproducible database use.
    """
    runtime = get_mash_runtime_version(mash_cmd)
    required = _version_tuple(db.mash_version)
    actual = _version_tuple(runtime)
    if actual[: len(required)] != required:
        raise PlasRipConfigError(
            f"Mash version mismatch: database requires {db.mash_version}, "
            f"but runtime '{mash_cmd}' reports {runtime}."
        )
    return runtime


def parse_mash_info_output(text: str) -> MashSketchInfo:
    """Parse the header block produced by ``mash info <file.msh>``."""
    patterns = {
        "hash_seed": re.compile(
            r"^\s*Hash function \(seed\):\s*(\S+)\s*\((\d+)\)\s*$", re.MULTILINE
        ),
        "k": re.compile(r"^\s*K-mer size:\s*(\d+)\b", re.MULTILINE),
        "alphabet": re.compile(r"^\s*Alphabet:\s*(.+?)\s*$", re.MULTILINE),
        "sketch_size": re.compile(
            r"^\s*Target min-hashes per sketch:\s*(\d+)\s*$", re.MULTILINE
        ),
        "sketches": re.compile(r"^\s*Sketches:\s*(\d+)\s*$", re.MULTILINE),
    }
    found = {name: pat.search(text) for name, pat in patterns.items()}
    missing = [name for name, match in found.items() if match is None]
    if missing:
        raise PlasRipConfigError(
            "Could not parse Mash sketch metadata from 'mash info' output; "
            f"missing fields: {missing}."
        )
    hash_function, seed = found["hash_seed"].groups()
    return MashSketchInfo(
        hash_function=hash_function,
        seed=int(seed),
        k=int(found["k"].group(1)),
        alphabet=found["alphabet"].group(1).strip(),
        sketch_size=int(found["sketch_size"].group(1)),
        sketches=int(found["sketches"].group(1)),
    )


def inspect_mash_sketch(
    msh_path: str | os.PathLike[str],
    mash_cmd: str = "mash",
) -> MashSketchInfo:
    """Inspect one Mash sketch using the runtime Mash executable."""
    path = Path(msh_path)
    if not path.is_file():
        raise PlasRipConfigError(f"Mash sketch not found: {path}")
    try:
        proc = subprocess.run(
            [mash_cmd, "info", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as exc:
        raise PlasRipConfigError(f"Mash executable not found: {mash_cmd!r}") from exc
    except subprocess.CalledProcessError as exc:
        detail = (exc.stderr or exc.stdout or "").strip()
        raise PlasRipConfigError(
            f"Failed to inspect Mash sketch {path}: {detail}"
        ) from exc
    return parse_mash_info_output((proc.stdout or "") + "\n" + (proc.stderr or ""))


def validate_mash_sketch_compatibility(
    db: PlasRipReferenceDatabase,
    info: MashSketchInfo,
    *,
    label: str = "Mash sketch",
) -> None:
    """Require a sketch to match all database-defining Mash parameters."""
    mismatches = []
    checks = [
        ("hash function", info.hash_function, db.mash_hash_function),
        ("seed", info.seed, db.mash_seed),
        ("k-mer size", info.k, db.mash_k),
        ("alphabet", info.alphabet, db.mash_alphabet),
        ("sketch size", info.sketch_size, db.mash_sketch_size),
        ("sketches per file", info.sketches, EXPECTED_SKETCHES_PER_FILE),
    ]
    for name, actual, expected in checks:
        if actual != expected:
            mismatches.append(f"{name}: expected {expected!r}, found {actual!r}")
    if mismatches:
        raise PlasRipConfigError(
            f"{label} is incompatible with database {db.database_version}: "
            + "; ".join(mismatches)
        )


def validate_reference_sketch_sample(
    db: PlasRipReferenceDatabase,
    mash_cmd: str = "mash",
    sample_size: int = 3,
) -> List[str]:
    """Inspect a deterministic sample of PTU representative sketches.

    This guards against a database directory containing sketches generated
    with parameters different from its manifest without paying the cost of
    invoking ``mash info`` for all 64,391 reference sketches on every run.
    """
    if sample_size < 1:
        return []
    index = db.plasmid_index.copy()
    index["_is_rep"] = _as_bool_series(index["is_PTU_representative"])
    reps = index.loc[index["_is_rep"], "plasmid_id"].astype(str).tolist()
    if not reps:
        raise PlasRipDataError("Reference database contains no PTU representative sketches.")

    if len(reps) <= sample_size:
        selected = reps
    else:
        positions = sorted(set(
            round(i * (len(reps) - 1) / (sample_size - 1))
            for i in range(sample_size)
        )) if sample_size > 1 else [0]
        selected = [reps[i] for i in positions]

    checked = []
    for plasmid_id in selected:
        path = db.mash_index_dir / f"{plasmid_id}.msh"
        info = inspect_mash_sketch(path, mash_cmd)
        validate_mash_sketch_compatibility(
            db, info, label=f"Reference sketch {path.name}"
        )
        checked.append(str(path))
    return checked

def build_assignment_indexes(
    db: PlasRipReferenceDatabase,
) -> Tuple[Dict[str, List[str]], List[str], Dict[str, Tuple[str, str]]]:
    """Build PTU member and representative indexes for assignment.

    Returns:
      ``ptu_to_members``: PTU -> stable member sketch paths
      ``rep_paths``: one representative sketch path per PTU
      ``rep_key_to_record``: Mash/reference-name variants -> (PTU_ID, plasmid_id)

    The PTU assignment algorithm itself remains implemented in ``assign.py``.
    """
    index = db.plasmid_index.copy()
    index["_is_rep"] = _as_bool_series(index["is_PTU_representative"])

    stable = index[
        (index["PTU_assignment_status"] == INDEX_ASSIGNED_STATUS)
        & (index["PTU_ID"] != "")
    ]

    ptu_to_members: Dict[str, List[str]] = {}
    missing_members: List[str] = []
    for row in stable.itertuples(index=False):
        ptu_id = str(getattr(row, "PTU_ID"))
        plasmid_id = str(getattr(row, "plasmid_id"))
        path = db.mash_index_dir / f"{plasmid_id}.msh"
        if path.is_file():
            ptu_to_members.setdefault(ptu_id, []).append(str(path))
        else:
            missing_members.append(str(path))

    if missing_members:
        raise PlasRipDataError(
            f"Reference database is incomplete: missing {len(missing_members)} PTU member sketches. "
            f"Examples: {missing_members[:5]}"
        )

    reps = index[index["_is_rep"]]
    rep_paths: List[str] = []
    rep_key_to_record: Dict[str, Tuple[str, str]] = {}
    for row in reps.itertuples(index=False):
        ptu_id = str(getattr(row, "PTU_ID"))
        plasmid_id = str(getattr(row, "plasmid_id"))
        mash_ref = str(getattr(row, "mash_reference_id"))
        path = db.mash_index_dir / f"{plasmid_id}.msh"
        if not path.is_file():
            raise PlasRipDataError(f"Representative Mash sketch not found: {path}")
        rep_paths.append(str(path))

        ref_base = Path(mash_ref).name
        ref_stem = Path(ref_base).stem
        keys = {
            mash_ref,
            ref_base,
            ref_stem,
            plasmid_id,
            f"{plasmid_id}.msh",
        }
        for key in keys:
            existing = rep_key_to_record.get(key)
            value = (ptu_id, plasmid_id)
            if existing is not None and existing != value:
                raise PlasRipDataError(
                    f"Ambiguous representative lookup key {key!r}: {existing} versus {value}"
                )
            rep_key_to_record[key] = value

    return ptu_to_members, rep_paths, rep_key_to_record


def ptu_metadata_lookup(db: PlasRipReferenceDatabase) -> Dict[str, Dict[str, str]]:
    """Return PTU metadata keyed by formal PTU ID."""
    return db.ptu_metadata.set_index("PTU_ID", drop=False).to_dict(orient="index")


def pgg_metadata_lookup(db: PlasRipReferenceDatabase) -> Dict[str, Dict[str, str]]:
    """Return PGG metadata keyed by formal PGG ID."""
    return db.pgg_metadata.set_index("PGG_ID", drop=False).to_dict(orient="index")
