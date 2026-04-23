"""Path helpers wrapping central config path functions."""

from __future__ import annotations

from pathlib import Path

from ..config import (
    get_conjscan_dir,
    get_data_dir,
    get_default_hostrange_dir,
    get_package_dir,
    get_ptu_info_path,
    get_ptu_member_path,
    get_ptu_metadata_path,
)

__all__ = [
    "get_package_dir",
    "get_data_dir",
    "get_default_hostrange_dir",
    "get_conjscan_dir",
    "get_ptu_metadata_path",
    "get_ptu_member_path",
    "get_ptu_info_path",
    "resolve_optional_path",
]


def resolve_optional_path(path: str | None) -> Path | None:
    """Resolve an optional path string to Path."""
    if path is None:
        return None
    return Path(path).resolve()
