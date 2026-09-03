"""Path helpers wrapping central config path functions."""

from __future__ import annotations

from pathlib import Path

from ..config import (
    get_conjscan_dir,
    get_data_dir,
    get_default_hostrange_dir,
    get_package_dir,
)

__all__ = [
    "get_package_dir",
    "get_data_dir",
    "get_default_hostrange_dir",
    "get_conjscan_dir",
    "resolve_optional_path",
]


def resolve_optional_path(path: str | None) -> Path | None:
    if path is None:
        return None
    return Path(path).resolve()
