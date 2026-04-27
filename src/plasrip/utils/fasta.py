"""FASTA validation helpers."""

from __future__ import annotations

from pathlib import Path

from ..exceptions import PlasRipInputError


_WILDCARD_CHARS = set("*?[")


def _looks_like_unexpanded_glob(path: str | Path) -> bool:
    return any(ch in str(path) for ch in _WILDCARD_CHARS)


def validate_fasta(path: str | Path) -> bool:
    """Return True if the path points to a readable FASTA with >=1 record."""
    p = Path(path).expanduser().resolve()
    if not p.is_file() or p.stat().st_size <= 0:
        return False
    try:
        with p.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith(">"):  # first record is enough for validation
                    return True
    except OSError:
        return False
    return False


def validate_fastas(paths: list[str | Path]) -> list[str]:
    """Validate multiple FASTA files and return absolute string paths."""
    if not paths:
        raise PlasRipInputError("No input FASTA file was provided.")

    valid: list[str] = []
    for path in paths:
        raw = str(path)
        p = Path(path).expanduser().resolve()

        if not p.exists():
            if _looks_like_unexpanded_glob(raw):
                raise PlasRipInputError(
                    f"Input FASTA path was not expanded by the shell: {raw}\n"
                    "No matching file was found. Quote handling or wildcard expansion may be wrong; "
                    "pass an existing FASTA file path or expand the wildcard before running PlasRiP."
                )
            raise PlasRipInputError(f"Input FASTA file not found: {p}")

        if p.is_dir():
            raise PlasRipInputError(f"Input FASTA path is a directory, expected a file: {p}")

        if p.stat().st_size <= 0:
            raise PlasRipInputError(f"Input FASTA file is empty: {p}")

        if not validate_fasta(p):
            raise PlasRipInputError(
                f"Invalid FASTA file: {p}\n"
                "No FASTA record header starting with '>' was found, or the file could not be read."
            )

        valid.append(str(p))
    return valid


def get_fasta_basename(path: str | Path) -> str:
    """Return the file name of a FASTA path."""
    return Path(path).name


def get_plasmid_id_from_fasta(path: str | Path) -> str:
    """Return the stem of a FASTA path."""
    return Path(path).stem
