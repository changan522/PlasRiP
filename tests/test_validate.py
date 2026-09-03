from __future__ import annotations

from pathlib import Path

import joblib
import pytest

from plasrip.exceptions import PlasRipConfigError, PlasRipModelError
from plasrip.utils.taxonomy import load_host_metadata
from plasrip.utils.validate import validate_model_dir, validate_out_dir


def test_validate_out_dir_creates_directory(tmp_path: Path):
    out_dir = tmp_path / "out"
    resolved = validate_out_dir(str(out_dir), overwrite=False)
    assert Path(resolved).is_dir()


def test_validate_out_dir_rejects_nonempty_directory_without_overwrite(tmp_path: Path):
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    (out_dir / "existing.txt").write_text("x", encoding="utf-8")
    with pytest.raises(PlasRipConfigError):
        validate_out_dir(str(out_dir), overwrite=False)


def test_validate_model_dir_structured_layout(tmp_path: Path):
    (tmp_path / "model").mkdir()
    (tmp_path / "cache").mkdir()
    (tmp_path / "host_db").mkdir()
    (tmp_path / "model" / "two_tower_best_3.2.pth").write_text("x", encoding="utf-8")
    (tmp_path / "cache" / "host_profiles.joblib").write_text("x", encoding="utf-8")
    (tmp_path / "cache" / "host_vocab.joblib").write_text("x", encoding="utf-8")
    (tmp_path / "cache" / "plasmid_codon_scaler.joblib").write_text("x", encoding="utf-8")
    (tmp_path / "host_db" / "host_metadata.tsv").write_text(
        "Host_id\tHost_species\n0\tsp1\n", encoding="utf-8"
    )

    report = validate_model_dir(str(tmp_path))

    assert Path(report["root"]).resolve() == tmp_path.resolve()
    assert Path(report["model_path"]).name == "two_tower_best_3.2.pth"
    assert report["model_name"] == "two_tower_best_3.2.pth"
    assert Path(report["host_profiles_path"]).name == "host_profiles.joblib"
    assert Path(report["host_vocab_path"]).name == "host_vocab.joblib"
    assert Path(report["plasmid_codon_scaler_path"]).name == "plasmid_codon_scaler.joblib"
    assert Path(report["host_metadata_path"]).name == "host_metadata.tsv"
    assert report["model/two_tower_best_3.2.pth"] is not None
    assert report["cache/host_profiles.joblib"] is not None
    assert report["cache/host_vocab.joblib"] is not None
    assert report["cache/plasmid_codon_scaler.joblib"] is not None


def test_validate_model_dir_flat_layout_and_taxonomy_fallback(tmp_path: Path):
    (tmp_path / "two_tower_best_3.2.pth").write_text("x", encoding="utf-8")
    (tmp_path / "host_profiles.joblib").write_text("x", encoding="utf-8")
    joblib.dump(
        {
            "host_key2id": {"escherichia_coli": 0, "klebsiella_pneumoniae": 1},
            "id2host_key": {0: "escherichia_coli", 1: "klebsiella_pneumoniae"},
        },
        tmp_path / "host_vocab.joblib",
    )
    (tmp_path / "plasmid_codon_scaler.joblib").write_text("x", encoding="utf-8")

    report = validate_model_dir(str(tmp_path))

    assert Path(report["root"]).resolve() == tmp_path.resolve()
    assert Path(report["model_path"]).name == "two_tower_best_3.2.pth"
    assert report["model_name"] == "two_tower_best_3.2.pth"
    assert Path(report["host_profiles_path"]).name == "host_profiles.joblib"
    assert Path(report["host_vocab_path"]).name == "host_vocab.joblib"
    assert Path(report["plasmid_codon_scaler_path"]).name == "plasmid_codon_scaler.joblib"
    assert report["host_metadata_path"] is None
    assert report["model/two_tower_best_3.2.pth"] is not None
    assert report["cache/host_profiles.joblib"] is not None
    assert report["cache/host_vocab.joblib"] is not None
    assert report["cache/plasmid_codon_scaler.joblib"] is not None

    meta = load_host_metadata(tmp_path)
    assert set(meta.columns) >= {"Host_id", "Host_species", "Host_genus", "Host_family"}
    assert list(meta["Host_species"]) == ["escherichia_coli", "klebsiella_pneumoniae"]
    assert list(meta["Host_genus"]) == ["escherichia", "klebsiella"]


def test_validate_model_dir_raises_on_missing_files(tmp_path: Path):
    with pytest.raises(PlasRipModelError):
        validate_model_dir(str(tmp_path))