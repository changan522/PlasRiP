"""Validation helpers for arguments, model dirs, and dependencies (3.2 unified)."""

from __future__ import annotations

import shutil
from pathlib import Path

from ..config import DEFAULT_MODEL_FILENAME, REQUIRED_MODEL_FILENAMES, get_hostrange_data_dir
from ..exceptions import PlasRipConfigError, PlasRipDependencyError, PlasRipModelError


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


def validate_model_dir(model_dir: str | None) -> dict[str, str | None]:
    if model_dir is None:
        model_dir = str(get_hostrange_data_dir())

    root = Path(model_dir).resolve()
    if not root.is_dir():
        raise PlasRipModelError(f"3.2 model directory not found: {root}")

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
            f"3.2 model directory is incomplete: {root}\n"
            f"Missing required files: {', '.join(missing)}"
        )

    return {
        "root": str(root),
        "model_path": model_path,
        "model_name": Path(model_path).name if model_path else None,
        "host_profiles_path": host_profiles_path,
        "host_vocab_path": host_vocab_path,
        "plasmid_codon_scaler_path": scaler_path,
        "host_metadata_path": host_meta_path,
        f"model/{DEFAULT_MODEL_FILENAME}": model_path,
        "cache/host_profiles.joblib": host_profiles_path,
        "cache/host_vocab.joblib": host_vocab_path,
        "cache/plasmid_codon_scaler.joblib": scaler_path,
    }


def validate_rep_dir(rep_dir: str) -> str:
    root = Path(rep_dir).resolve()
    if not root.is_dir():
        raise PlasRipConfigError(f"PTU sketch directory not found: {root}")
    return str(root)


def validate_out_dir(out_dir: str, overwrite: bool = False) -> str:
    root = Path(out_dir).resolve()
    if root.exists() and any(root.iterdir()) and not overwrite:
        raise PlasRipConfigError(
            f"Output directory already exists and is not empty: {root}. Use overwrite=True to reuse it."
        )
    root.mkdir(parents=True, exist_ok=True)
    return str(root)
