"""FASTA validation helpers."""

from __future__ import annotations

from pathlib import Path

from ..exceptions import PlasRipInputError


def validate_fasta(path: str | Path) -> bool:
    """Return True if the path points to a readable FASTA with >=1 record."""
    p = Path(path).resolve()
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
    valid: list[str] = []
    for path in paths:
        p = Path(path).resolve()
        if not validate_fasta(p):
            raise PlasRipInputError(f"Invalid FASTA file: {p}")
        valid.append(str(p))
    return valid


def get_fasta_basename(path: str | Path) -> str:
    """Return the file name of a FASTA path."""
    return Path(path).name


def get_plasmid_id_from_fasta(path: str | Path) -> str:
    """Return the stem of a FASTA path."""
    return Path(path).stem
