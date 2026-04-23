from __future__ import annotations

from pathlib import Path

from plasrip.config import get_data_dir, get_package_dir
from plasrip.utils.io import ensure_dir, write_tsv


def test_package_paths_exist():
    assert get_package_dir().is_dir()
    assert get_data_dir().is_dir() or not get_data_dir().exists()  # package may be prepared before data are copied


def test_write_tsv_smoke(tmp_path: Path):
    out_dir = ensure_dir(tmp_path / "outputs", overwrite=True)
    out_path = write_tsv(__import__("pandas").DataFrame([{"a": 1}]), out_dir / "test.tsv")
    assert out_path.is_file()
