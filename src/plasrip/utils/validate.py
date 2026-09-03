"""Validation helpers for arguments, model dirs, and dependencies (3.2 unified)."""

from __future__ import annotations

import shutil
from pathlib import Path

from ..config import (
    DEFAULT_MODEL_FILENAME,
    REQUIRED_MODEL_FILENAMES,
    get_env_hostrange_dir,
    get_hostrange_data_dir,
)
from ..exceptions import PlasRipConfigError, PlasRipDependencyError, PlasRipModelError
from ..database import validate_database_dir as _validate_database_dir


def check_external_dependency(cmd: str) -> bool:
    return shutil.which(cmd) is not None


def check_dependencies(
    need_mash: bool = False,
    need_abricate: bool = False,
    need_prodigal: bool = False,
    need_macsyfinder: bool = False,
) -> dict[str, bool]:
    checks = {
        "mash": (need_mash, check_external_dependency("mash")),
        "abricate": (need_abricate, check_external_dependency("abricate")),
        "prodigal": (need_prodigal, check_external_dependency("prodigal")),
        "macsyfinder": (need_macsyfinder, check_external_dependency("macsyfinder")),
    }
    for name, (needed, ok) in checks.items():
        if needed and not ok:
            raise PlasRipDependencyError(f"Required dependency not found on PATH: {name}")
    return {name: ok for name, (_, ok) in checks.items()}


def _pick_existing(root: Path, rels: list[Path]) -> str | None:
    for rel in rels:
        p = root / rel
        if p.is_file():
            return str(p.resolve())
    return None


def _format_model_dir_hint(root: Path | None = None) -> str:
    bundled = get_hostrange_data_dir()
    parts = [
        "Host-range model directory could not be validated.",
        f"Bundled default location: {bundled}",
        "Expected files: "
        f"{DEFAULT_MODEL_FILENAME}, host_profiles.joblib, host_vocab.joblib, "
        "plasmid_codon_scaler.joblib",
    ]
    if root is not None:
        parts.insert(1, f"Checked location: {root}")
    parts.append("Use --model-dir <DIR> or set PLASRIP_MODEL_DIR=<DIR> if your models are stored elsewhere.")
    return "\n".join(parts)


def resolve_model_dir(model_dir: str | None = None) -> Path:
    """Resolve the host-range model directory.

    Priority:
    1. explicit --model-dir
    2. PLASRIP_MODEL_DIR
    3. bundled package data directory, e.g.
       <env>/lib/pythonX/site-packages/plasrip/data/hostrange
       or <source>/src/plasrip/data/hostrange for editable/source runs
    """
    if model_dir:
        return Path(model_dir).expanduser().resolve()

    env_dir = get_env_hostrange_dir()
    if env_dir is not None:
        return env_dir

    return get_hostrange_data_dir().resolve()


def validate_model_dir(model_dir: str | None) -> dict[str, str | None]:
    root = resolve_model_dir(model_dir)

    if not root.exists():
        raise PlasRipModelError(
            f"Host-range model directory not found: {root}\n"
            + _format_model_dir_hint(root)
        )
    if not root.is_dir():
        raise PlasRipModelError(
            f"Host-range model path is not a directory: {root}\n"
            + _format_model_dir_hint(root)
        )

    model_candidates = [Path("model") / DEFAULT_MODEL_FILENAME, Path(DEFAULT_MODEL_FILENAME)]
    model_path = _pick_existing(root, model_candidates)
    host_profiles_path = _pick_existing(root, [Path("cache") / "host_profiles.joblib", Path("host_profiles.joblib")])
    host_vocab_path = _pick_existing(root, [Path("cache") / "host_vocab.joblib", Path("host_vocab.joblib")])
    scaler_path = _pick_existing(root, [Path("cache") / "plasmid_codon_scaler.joblib", Path("plasmid_codon_scaler.joblib")])
    host_meta_path = _pick_existing(root, [Path("host_db") / "host_metadata.tsv", Path("host_metadata.tsv")])

    missing: list[str] = []
    if model_path is None:
        missing.extend(REQUIRED_MODEL_FILENAMES)
    if host_profiles_path is None:
        missing.append("host_profiles.joblib")
    if host_vocab_path is None:
        missing.append("host_vocab.joblib")
    if scaler_path is None:
        missing.append("plasmid_codon_scaler.joblib")

    if missing:
        raise PlasRipModelError(
            f"Host-range model directory is incomplete: {root}\n"
            f"Missing required files: {', '.join(missing)}\n"
            + _format_model_dir_hint(root)
        )

    return {
        "root": str(root),
        "model_path": model_path,
        "model_name": Path(model_path).name if model_path else None,
        "host_profiles_path": host_profiles_path,
        "host_vocab_path": host_vocab_path,
        "plasmid_codon_scaler_path": scaler_path,
        "host_metadata_path": host_meta_path,
        # Backward-compatible structured keys used by older reports/tests.
        f"model/{DEFAULT_MODEL_FILENAME}": model_path,
        "cache/host_profiles.joblib": host_profiles_path,
        "cache/host_vocab.joblib": host_vocab_path,
        "cache/plasmid_codon_scaler.joblib": scaler_path,
    }

def validate_database_dir(database_dir: str | None, check_sketches: bool = False) -> dict[str, object]:
    """Validate the external PlasRiP reference database."""
    return _validate_database_dir(database_dir, check_sketches=check_sketches)


def validate_out_dir(out_dir: str, overwrite: bool = False) -> str:
    root = Path(out_dir).resolve()
    if root.exists() and any(root.iterdir()) and not overwrite:
        raise PlasRipConfigError(
            f"Output directory already exists and is not empty: {root}. Use overwrite=True to reuse it."
        )
    root.mkdir(parents=True, exist_ok=True)
    return str(root)
