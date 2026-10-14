from __future__ import annotations

from pathlib import Path
import pandas as pd
import pytest

from plasrip.assign import parse_rep_hits
from plasrip.database import (
    INDEX_REQUIRED_COLUMNS,
    PGG_REQUIRED_COLUMNS,
    MPC_REQUIRED_COLUMNS,
    build_assignment_indexes,
    load_database,
    resolve_database_dir,
    validate_database_dir,
    validate_query_mash_parameters,
)
from plasrip.exceptions import PlasRipConfigError, PlasRipDataError


def _complete_row(columns: list[str], values: dict[str, str]) -> dict[str, str]:
    row = {c: "" for c in columns}
    row.update(values)
    return row


def make_database(root: Path) -> Path:
    (root / "mash_index").mkdir(parents=True)

    index_rows = [
        _complete_row(INDEX_REQUIRED_COLUMNS, {
            "mash_reference_id": "plasmid_fasta/ref1.fasta",
            "plasmid_id": "ref1",
            "accession": "A1",
            "length_bp": "1000",
            "MPC_ID": "MPC000001",
            "MPC_assignment_status": "ASSIGNED_STABLE_MPC",
            "is_MPC_representative": "TRUE",
        }),
        _complete_row(INDEX_REQUIRED_COLUMNS, {
            "mash_reference_id": "plasmid_fasta/ref2.fasta",
            "plasmid_id": "ref2",
            "accession": "A2",
            "length_bp": "900",
            "MPC_ID": "MPC000001",
            "MPC_assignment_status": "ASSIGNED_STABLE_MPC",
            "is_MPC_representative": "FALSE",
        }),
    ]
    pd.DataFrame(index_rows).to_csv(root / "plasmid_reference_index.tsv", sep="\t", index=False)

    mpc_row = _complete_row(MPC_REQUIRED_COLUMNS, {
        "MPC_ID": "MPC000001",
        "representative_plasmid": "ref1",
        "n_plasmids": "2",
        "PGG_ID": "PGG-000001",
        "PGG_assignment_status": "ASSIGNED_STABLE_PGG",
        "mean_length_bp": "950",
        "median_length_bp": "950",
        "min_length_bp": "900",
        "max_length_bp": "1000",
        "replicon_annotated_fraction": "1",
        "n_replicon_types": "1",
        "replicon_types": "IncX",
        "dominant_replicon": "IncX",
        "dominant_replicon_fraction": "1",
        "n_host_species": "1",
        "observed_host_species": "Species A",
        "top_host_species": "Species A",
        "top_host_species_count": "2",
        "top_host_species_coverage": "1",
        "n_host_genera": "1",
        "observed_host_genera": "GenusA",
        "top_host_genus": "GenusA",
        "top_host_genus_count": "2",
        "top_host_genus_coverage": "1",
        "n_host_families": "1",
        "top_host_family": "FamilyA",
        "top_host_family_coverage": "1",
        "n_host_orders": "1",
        "top_host_order": "OrderA",
        "top_host_order_coverage": "1",
        "n_host_classes": "1",
        "top_host_class": "ClassA",
        "top_host_class_coverage": "1",
        "n_host_phyla": "1",
        "top_host_phylum": "PhylumA",
        "top_host_phylum_count": "2",
        "top_host_phylum_coverage": "1",
        "observed_host_range_level": "1",
        "observed_host_range_rank": "Species",
        "dominant_taxon": "Species A",
        "dominant_taxon_coverage": "1",
        "n_CONJ": "1",
        "CONJ_fraction": "0.5",
        "n_MOB": "1",
        "MOB_fraction": "0.5",
        "n_NT": "0",
        "NT_fraction": "0",
        "dominant_mobility": "CONJ",
        "n_ARG_positive_plasmids": "1",
        "ARG_positive_fraction": "0.5",
        "mean_ARG_count_per_plasmid": "0.5",
        "n_unique_ARGs": "1",
        "unique_ARGs": "blaX",
        "n_unique_ARG_classes": "1",
        "unique_ARG_classes": "Beta-lactam",
        "n_virulence_positive_plasmids": "1",
        "virulence_positive_fraction": "0.5",
        "mean_virulence_gene_count_per_plasmid": "0.5",
        "n_unique_virulence_genes": "1",
        "unique_virulence_genes": "virX",
    })
    pd.DataFrame([mpc_row]).to_csv(root / "MPC_metadata.tsv", sep="\t", index=False)

    pgg_row = _complete_row(PGG_REQUIRED_COLUMNS, {
        "PGG_ID": "PGG-000001",
        "n_MPCs": "1",
        "n_plasmids": "2",
        "anchor_MPC": "MPC000001",
        "anchor_representative_plasmid": "ref1",
        "minimum_internal_support": "0.9",
        "merge_with_parent_support": "0.5",
        "effective_persistence_width": "0.4",
        "persistence_volume": "0.4",
        "anchor_mean_within_group_Mash": "0.01",
    })
    pd.DataFrame([pgg_row]).to_csv(root / "PGG_metadata.tsv", sep="\t", index=False)

    pd.DataFrame([
        ("database_name", "PlasRiP Reference Database"),
        ("database_version", "v3.0-test"),
        ("schema_version", "2.0"),
        ("release_date", "2026-09-03"),
        ("mash_index_directory", "mash_index/"),
        ("total_reference_plasmids", "2"),
        ("stable_MPC_assigned_plasmids", "2"),
        ("unassigned_no_stable_MPC_plasmids", "0"),
        ("MPC_count", "1"),
        ("PGG_assigned_MPCs", "1"),
        ("PGG_count", "1"),
        ("MPC_ID_format", "MPC000001"),
        ("PGG_ID_format", "PGG-000001"),
        ("mash_version", "2.3"),
        ("mash_hash_function", "MurmurHash3_x64_128"),
        ("mash_seed", "42"),
        ("mash_k", "21"),
        ("mash_alphabet", "ACGT (canonical)"),
        ("mash_target_min_hashes_per_sketch", "1000"),
        ("mash_sketches_per_file", "1"),
        ("mash_reference_id_pattern", "plasmid_fasta/<sequence_filename>"),
        ("query_MPC_Mash_cutoff", "0.05"),
    ], columns=["key", "value"]).to_csv(root / "database_manifest.tsv", sep="\t", index=False)

    for name in ["ref1.msh", "ref2.msh"]:
        (root / "mash_index" / name).write_bytes(b"test")
    return root


def _rewrite_tsv(path: Path, mutator):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    mutator(df)
    df.to_csv(path, sep="\t", index=False)


def _set_manifest_value(path: Path, key: str, value: str) -> None:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df.loc[df["key"] == key, "value"] = value
    df.to_csv(path, sep="\t", index=False)


def test_database_load_and_indexes(tmp_path: Path):
    db_dir = make_database(tmp_path / "db")
    report = validate_database_dir(db_dir, check_sketches=True, check_all_sketches=True)
    assert report["MPCs"] == 1
    assert report["PGGs"] == 1
    assert report["mash_k"] == 21
    assert report["mash_sketch_size"] == 1000
    assert report["mash_alphabet"] == "ACGT (canonical)"

    db = load_database(db_dir)
    assert db.schema_version == "2.0"
    assert db.mash_version == "2.3"
    assert db.mash_seed == 42

    mpc_to_members, reps, keys = build_assignment_indexes(db)
    assert len(mpc_to_members["MPC000001"]) == 2
    assert len(reps) == 1
    assert keys["plasmid_fasta/ref1.fasta"] == ("MPC000001", "ref1")


def test_resolve_database_from_environment(tmp_path: Path, monkeypatch):
    db_dir = make_database(tmp_path / "db")
    monkeypatch.setenv("PLASRIP_DATABASE", str(db_dir))
    assert resolve_database_dir(None) == db_dir.resolve()
    assert validate_database_dir(None)["MPCs"] == 1


def test_explicit_database_takes_priority_over_environment(tmp_path: Path, monkeypatch):
    db1 = make_database(tmp_path / "db1")
    db2 = make_database(tmp_path / "db2")
    monkeypatch.setenv("PLASRIP_DATABASE", str(db2))
    assert resolve_database_dir(db1) == db1.resolve()


def test_parse_rep_hits_uses_mash_reference_id():
    rec = {"plasmid_fasta/ref1.fasta": ("MPC000001", "ref1")}
    lines = ["query.fasta\tplasmid_fasta/ref1.fasta\t0.01\t0\t1/1"]
    hits = parse_rep_hits("query.fasta", lines, rec, 0.05)
    assert hits[0]["MPC_ID"] == "MPC000001"
    assert hits[0]["Representative"] == "ref1"


def test_query_mash_parameters_must_match_database(tmp_path: Path):
    db = load_database(make_database(tmp_path / "db"))
    assert validate_query_mash_parameters(db, None, None) == (21, 1000)
    with pytest.raises(PlasRipConfigError):
        validate_query_mash_parameters(db, 31, None)
    with pytest.raises(PlasRipConfigError):
        validate_query_mash_parameters(db, None, 5000)


def test_database_rejects_member_count_mismatch(tmp_path: Path):
    db_dir = make_database(tmp_path / "db")
    _rewrite_tsv(db_dir / "MPC_metadata.tsv", lambda df: df.__setitem__("n_plasmids", "99"))
    with pytest.raises(PlasRipDataError, match="member-count mismatch"):
        validate_database_dir(db_dir)


def test_database_rejects_representative_mismatch(tmp_path: Path):
    db_dir = make_database(tmp_path / "db")
    _rewrite_tsv(
        db_dir / "MPC_metadata.tsv",
        lambda df: df.__setitem__("representative_plasmid", "ref2"),
    )
    with pytest.raises(PlasRipDataError, match="representative mismatch"):
        validate_database_dir(db_dir)


def test_database_rejects_unsupported_schema(tmp_path: Path):
    db_dir = make_database(tmp_path / "db")
    _set_manifest_value(db_dir / "database_manifest.tsv", "schema_version", "9.9")
    with pytest.raises(PlasRipDataError, match="Unsupported database schema_version"):
        validate_database_dir(db_dir)


def test_database_rejects_missing_schema_column(tmp_path: Path):
    db_dir = make_database(tmp_path / "db")
    df = pd.read_csv(db_dir / "MPC_metadata.tsv", sep="\t", dtype=str, keep_default_na=False)
    df = df.drop(columns=["unique_ARG_classes"])
    df.to_csv(db_dir / "MPC_metadata.tsv", sep="\t", index=False)
    with pytest.raises(PlasRipDataError, match="missing required schema-2.0 columns"):
        validate_database_dir(db_dir)


def test_database_rejects_invalid_boolean_encoding(tmp_path: Path):
    db_dir = make_database(tmp_path / "db")
    _rewrite_tsv(
        db_dir / "plasmid_reference_index.tsv",
        lambda df: df.loc[df["plasmid_id"] == "ref2", "is_MPC_representative"].__setitem__(slice(None), "NOPE"),
    )
    # Avoid chained-assignment ambiguity across pandas versions.
    df = pd.read_csv(db_dir / "plasmid_reference_index.tsv", sep="\t", dtype=str, keep_default_na=False)
    df.loc[df["plasmid_id"] == "ref2", "is_MPC_representative"] = "NOPE"
    df.to_csv(db_dir / "plasmid_reference_index.tsv", sep="\t", index=False)
    with pytest.raises(PlasRipDataError, match="canonical TRUE/FALSE"):
        validate_database_dir(db_dir)


def test_database_rejects_mash_index_escape(tmp_path: Path):
    db_dir = make_database(tmp_path / "db")
    _set_manifest_value(db_dir / "database_manifest.tsv", "mash_index_directory", "../outside")
    with pytest.raises(PlasRipDataError, match="escapes the database root"):
        validate_database_dir(db_dir)


def test_database_rejects_bad_mash_manifest_semantics(tmp_path: Path):
    db_dir = make_database(tmp_path / "db")
    _set_manifest_value(db_dir / "database_manifest.tsv", "mash_sketches_per_file", "2")
    with pytest.raises(PlasRipDataError, match="exactly 1 sketch"):
        validate_database_dir(db_dir)


def test_database_rejects_host_range_level_rank_mismatch(tmp_path: Path):
    db_dir = make_database(tmp_path / "db")
    _rewrite_tsv(
        db_dir / "MPC_metadata.tsv",
        lambda df: df.__setitem__("observed_host_range_rank", "Genus"),
    )
    with pytest.raises(PlasRipDataError, match="host_range_level/rank mismatch"):
        validate_database_dir(db_dir)


def test_database_rejects_mobility_count_mismatch(tmp_path: Path):
    db_dir = make_database(tmp_path / "db")
    _rewrite_tsv(db_dir / "MPC_metadata.tsv", lambda df: df.__setitem__("n_NT", "1"))
    with pytest.raises(PlasRipDataError, match="mobility counts do not sum"):
        validate_database_dir(db_dir)


def test_assignment_indexes_exclude_unassigned_reference_plasmids(tmp_path: Path):
    db_dir = make_database(tmp_path / "db")

    index = pd.read_csv(
        db_dir / "plasmid_reference_index.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    index.loc[len(index)] = {
        "mash_reference_id": "plasmid_fasta/unassigned.fasta",
        "plasmid_id": "unassigned",
        "accession": "A3",
        "length_bp": "800",
        "MPC_ID": "",
        "MPC_assignment_status": "UNASSIGNED_NO_STABLE_MPC",
        "is_MPC_representative": "FALSE",
    }
    index.to_csv(db_dir / "plasmid_reference_index.tsv", sep="\t", index=False)

    manifest = pd.read_csv(
        db_dir / "database_manifest.tsv", sep="\t", dtype=str, keep_default_na=False
    )
    manifest.loc[manifest["key"] == "total_reference_plasmids", "value"] = "3"
    manifest.loc[manifest["key"] == "unassigned_no_stable_MPC_plasmids", "value"] = "1"
    manifest.to_csv(db_dir / "database_manifest.tsv", sep="\t", index=False)

    # Normal assignment only requires stable MPC member/repr sketches; the
    # unassigned reference remains indexed but is not a MPC-assignment target.
    db = load_database(db_dir)
    mpc_to_members, reps, keys = build_assignment_indexes(db)
    assert len(mpc_to_members["MPC000001"]) == 2
    assert len(reps) == 1
    assert "unassigned" not in keys


def test_parse_mash_info_output_matches_database_header():
    from plasrip.database import parse_mash_info_output

    text = """Header:\n  Hash function (seed):          MurmurHash3_x64_128 (42)\n  K-mer size:                    21 (64-bit hashes)\n  Alphabet:                      ACGT (canonical)\n  Target min-hashes per sketch:  1000\n  Sketches:                      1\n\nSketches:\n"""
    info = parse_mash_info_output(text)
    assert info.hash_function == "MurmurHash3_x64_128"
    assert info.seed == 42
    assert info.k == 21
    assert info.alphabet == "ACGT (canonical)"
    assert info.sketch_size == 1000
    assert info.sketches == 1


def test_mash_runtime_version_must_match_database_major_minor(tmp_path: Path, monkeypatch):
    from plasrip import database as dbmod

    db = load_database(make_database(tmp_path / "db"))
    monkeypatch.setattr(dbmod, "get_mash_runtime_version", lambda *a, **k: "2.3")
    assert dbmod.validate_mash_runtime(db, "mash") == "2.3"

    monkeypatch.setattr(dbmod, "get_mash_runtime_version", lambda *a, **k: "2.3.1")
    assert dbmod.validate_mash_runtime(db, "mash") == "2.3.1"

    monkeypatch.setattr(dbmod, "get_mash_runtime_version", lambda *a, **k: "2.2")
    with pytest.raises(PlasRipConfigError, match="version mismatch"):
        dbmod.validate_mash_runtime(db, "mash")


def test_mash_sketch_compatibility_checks_seed_k_size_and_alphabet(tmp_path: Path):
    from plasrip.database import MashSketchInfo, validate_mash_sketch_compatibility

    db = load_database(make_database(tmp_path / "db"))
    good = MashSketchInfo(
        hash_function="MurmurHash3_x64_128",
        seed=42,
        k=21,
        alphabet="ACGT (canonical)",
        sketch_size=1000,
        sketches=1,
    )
    validate_mash_sketch_compatibility(db, good)

    bad = MashSketchInfo(
        hash_function="MurmurHash3_x64_128",
        seed=7,
        k=21,
        alphabet="ACGT (canonical)",
        sketch_size=1000,
        sketches=1,
    )
    with pytest.raises(PlasRipConfigError, match="seed"):
        validate_mash_sketch_compatibility(db, bad, label="query")


def test_reference_sketch_sample_is_deterministic(tmp_path: Path, monkeypatch):
    from plasrip import database as dbmod
    from plasrip.database import MashSketchInfo

    db_dir = make_database(tmp_path / "db")
    db = load_database(db_dir)
    seen = []

    def fake_inspect(path, mash_cmd):
        seen.append(Path(path).name)
        return MashSketchInfo(
            hash_function="MurmurHash3_x64_128",
            seed=42,
            k=21,
            alphabet="ACGT (canonical)",
            sketch_size=1000,
            sketches=1,
        )

    monkeypatch.setattr(dbmod, "inspect_mash_sketch", fake_inspect)
    checked = dbmod.validate_reference_sketch_sample(db, "mash", sample_size=3)
    assert len(checked) == 1  # one MPC representative in this fixture
    assert seen == ["ref1.msh"]
