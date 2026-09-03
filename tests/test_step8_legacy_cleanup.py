from __future__ import annotations

from pathlib import Path
import tomllib


LEGACY_PTU_FILES = {
    "PTU_info.tsv",
    "PTUmember.csv",
    "PTU-metadata.csv",
}


def _project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def test_legacy_ptu_database_files_are_not_bundled():
    data_dir = _project_root() / "src" / "plasrip" / "data"
    found = {p.name for p in data_dir.rglob("*") if p.is_file()}
    assert LEGACY_PTU_FILES.isdisjoint(found)


def test_package_data_is_explicit_and_excludes_external_reference_database():
    pyproject = tomllib.loads((_project_root() / "pyproject.toml").read_text(encoding="utf-8"))
    setuptools_cfg = pyproject["tool"]["setuptools"]
    assert setuptools_cfg["include-package-data"] is False

    patterns = set(setuptools_cfg["package-data"]["plasrip"])
    assert patterns == {"data/CONJScan/**", "data/hostrange/**"}


def test_runtime_source_does_not_name_legacy_ptu_database_files():
    src_dir = _project_root() / "src" / "plasrip"
    legacy_tokens = tuple(LEGACY_PTU_FILES)
    offenders = []
    for path in src_dir.rglob("*.py"):
        text = path.read_text(encoding="utf-8")
        if any(token in text for token in legacy_tokens):
            offenders.append(path.relative_to(_project_root()).as_posix())
    assert offenders == []
