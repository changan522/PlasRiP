"""General I/O helpers for PlasRiP v2."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


def ensure_dir(path: str | Path, overwrite: bool = False) -> Path:
    """Ensure a directory exists and optionally allow reusing a non-empty path."""
    p = Path(path).resolve()
    if p.exists() and not p.is_dir():
        raise NotADirectoryError(str(p))
    if p.exists() and any(p.iterdir()) and not overwrite:
        raise FileExistsError(
            f"Directory already exists and is not empty: {p}. "
            "Use overwrite=True to reuse it."
        )
    p.mkdir(parents=True, exist_ok=True)
    return p


def ensure_parent_dir(path: str | Path) -> Path:
    """Create the parent directory for a file path and return it."""
    p = Path(path).resolve()
    p.parent.mkdir(parents=True, exist_ok=True)
    return p.parent


def ensure_file_exists(path: str | Path) -> Path:
    """Raise FileNotFoundError if the file does not exist."""
    p = Path(path).resolve()
    if not p.is_file():
        raise FileNotFoundError(str(p))
    return p


def write_tsv(df: pd.DataFrame, path: str | Path) -> Path:
    """Write a DataFrame as TSV and return the resolved path."""
    p = Path(path).resolve()
    ensure_parent_dir(p)
    df.to_csv(p, sep="\t", index=False)
    return p


def append_tsv(df: pd.DataFrame, path: str | Path) -> Path:
    """Append a DataFrame to a TSV file, writing the header only once."""
    p = Path(path).resolve()
    ensure_parent_dir(p)
    if df is None or df.empty:
        return p
    header = not p.exists() or p.stat().st_size == 0
    df.to_csv(p, sep="\t", index=False, mode="a", header=header)
    return p


def reset_tsv(path: str | Path) -> Path:
    """Remove a TSV output file if it exists and return its resolved path."""
    p = Path(path).resolve()
    if p.exists():
        p.unlink()
    return p


def read_tsv(path: str | Path) -> pd.DataFrame:
    """Read a TSV file into a DataFrame."""
    p = ensure_file_exists(path)
    return pd.read_csv(p, sep="\t")


def safe_remove(path: str | Path) -> None:
    """Remove a file if it exists."""
    p = Path(path).resolve()
    if p.exists():
        p.unlink()
