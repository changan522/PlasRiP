from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

import pandas as pd

from plasrip import annotate
from plasrip.schemas import ANNOTATION_SUMMARY_COLUMNS
from plasrip.database import INDEX_REQUIRED_COLUMNS, PTU_REQUIRED_COLUMNS, PGG_REQUIRED_COLUMNS


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
            "PTU_ID": "PTU000001",
            "PTU_assignment_status": "ASSIGNED_STABLE_PTU",
            "is_PTU_representative": "TRUE",
        }),
        _complete_row(INDEX_REQUIRED_COLUMNS, {
            "mash_reference_id": "plasmid_fasta/ref2.fasta",
            "plasmid_id": "ref2",
            "accession": "A2",
            "length_bp": "900",
            "PTU_ID": "PTU000001",
            "PTU_assignment_status": "ASSIGNED_STABLE_PTU",
            "is_PTU_representative": "FALSE",
        }),
    ]
    pd.DataFrame(index_rows).to_csv(root / "plasmid_reference_index.tsv", sep="\t", index=False)
    ptu_row = _complete_row(PTU_REQUIRED_COLUMNS, {
        "PTU_ID": "PTU000001", "representative_plasmid": "ref1", "n_plasmids": "2",
        "PGG_ID": "PGG-000001", "PGG_assignment_status": "ASSIGNED_STABLE_PGG",
        "mean_length_bp": "950", "median_length_bp": "950", "min_length_bp": "900", "max_length_bp": "1000",
        "replicon_annotated_fraction": "1", "n_replicon_types": "1", "replicon_types": "IncX",
        "dominant_replicon": "IncX", "dominant_replicon_fraction": "1",
        "n_host_species": "1", "observed_host_species": "Species A", "top_host_species": "Species A",
        "top_host_species_count": "2", "top_host_species_coverage": "1",
        "n_host_genera": "1", "observed_host_genera": "GenusA", "top_host_genus": "GenusA",
        "top_host_genus_count": "2", "top_host_genus_coverage": "1",
        "n_host_families": "1", "top_host_family": "FamilyA", "top_host_family_coverage": "1",
        "n_host_orders": "1", "top_host_order": "OrderA", "top_host_order_coverage": "1",
        "n_host_classes": "1", "top_host_class": "ClassA", "top_host_class_coverage": "1",
        "n_host_phyla": "1", "top_host_phylum": "PhylumA", "top_host_phylum_count": "2", "top_host_phylum_coverage": "1",
        "observed_host_range_level": "1", "observed_host_range_rank": "Species",
        "dominant_taxon": "Species A", "dominant_taxon_coverage": "1",
        "n_CONJ": "1", "CONJ_fraction": "0.5", "n_MOB": "1", "MOB_fraction": "0.5", "n_NT": "0", "NT_fraction": "0",
        "dominant_mobility": "CONJ",
        "n_ARG_positive_plasmids": "1", "ARG_positive_fraction": "0.5", "mean_ARG_count_per_plasmid": "0.5",
        "n_unique_ARGs": "1", "unique_ARGs": "blaX", "n_unique_ARG_classes": "1", "unique_ARG_classes": "Beta-lactam",
        "n_virulence_positive_plasmids": "1", "virulence_positive_fraction": "0.5",
        "mean_virulence_gene_count_per_plasmid": "0.5", "n_unique_virulence_genes": "1", "unique_virulence_genes": "virX",
    })
    pd.DataFrame([ptu_row]).to_csv(root / "PTU_metadata.tsv", sep="\t", index=False)
    pgg_row = _complete_row(PGG_REQUIRED_COLUMNS, {
        "PGG_ID": "PGG-000001", "n_PTUs": "1", "n_plasmids": "2", "anchor_PTU": "PTU000001",
        "anchor_representative_plasmid": "ref1", "minimum_internal_support": "0.9",
        "merge_with_parent_support": "0.5", "effective_persistence_width": "0.4",
        "persistence_volume": "0.4", "anchor_mean_within_group_Mash": "0.01",
    })
    pd.DataFrame([pgg_row]).to_csv(root / "PGG_metadata.tsv", sep="\t", index=False)
    pd.DataFrame([
        ("database_name", "PlasRiP Reference Database"), ("database_version", "v2.0-test"),
        ("schema_version", "2.0"), ("release_date", "2026-09-03"), ("mash_index_directory", "mash_index/"),
        ("total_reference_plasmids", "2"), ("stable_PTU_assigned_plasmids", "2"),
        ("unassigned_no_stable_PTU_plasmids", "0"), ("PTU_count", "1"), ("PGG_assigned_PTUs", "1"), ("PGG_count", "1"),
        ("PTU_ID_format", "PTU000001"), ("PGG_ID_format", "PGG-000001"), ("mash_version", "2.3"),
        ("mash_hash_function", "MurmurHash3_x64_128"), ("mash_seed", "42"), ("mash_k", "21"),
        ("mash_alphabet", "ACGT (canonical)"), ("mash_target_min_hashes_per_sketch", "1000"),
        ("mash_sketches_per_file", "1"), ("mash_reference_id_pattern", "plasmid_fasta/<sequence_filename>"),
        ("query_PTU_Mash_cutoff", "0.05"),
    ], columns=["key", "value"]).to_csv(root / "database_manifest.tsv", sep="\t", index=False)
    for name in ["ref1.msh", "ref2.msh"]:
        (root / "mash_index" / name).write_bytes(b"test")
    return root


def _write_fake_mash(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    script = path / "mash"
    script.write_text(
        r'''#!/usr/bin/env python3
import sys
from pathlib import Path

args = sys.argv[1:]
if args == ["--version"]:
    print("Mash version 2.3")
    raise SystemExit(0)

if args and args[0] == "info":
    print("""Header:
  Hash function (seed):          MurmurHash3_x64_128 (42)
  K-mer size:                    21 (64-bit hashes)
  Alphabet:                      ACGT (canonical)
  Target min-hashes per sketch:  1000
  Sketches:                      1

Sketches:
  [Hashes]  [Length]  [ID]  [Comment]
  1000      1000      plasmid_fasta/fake.fasta  fake
""")
    raise SystemExit(0)

if args and args[0] == "sketch":
    prefix = args[args.index("-o") + 1]
    Path(prefix + ".msh").write_bytes(b"fake-query-sketch")
    raise SystemExit(0)

if args and args[0] == "dist":
    xs = args[1:]
    if xs and xs[0] == "-p":
        xs = xs[2:]
    query = xs[0]
    refs = xs[1:]
    for ref in refs:
        base = Path(ref).stem
        dist = 0.01 if base == "ref1" else 0.02 if base == "ref2" else 0.20
        print(f"{query}\tplasmid_fasta/{base}.fasta\t{dist}\t0\t1000/1000")
    raise SystemExit(0)

print(f"unsupported fake Mash invocation: {args}", file=sys.stderr)
raise SystemExit(2)
''',
        encoding="utf-8",
    )
    script.chmod(0o755)
    return script


def test_step9_public_cli_database_and_assignment_subprocess(tmp_path: Path):
    """Exercise public argparse -> database -> Mash -> PTU/PGG TSV boundaries."""
    repo_root = Path(__file__).resolve().parents[1]
    db = make_database(tmp_path / "db")
    fake_mash = _write_fake_mash(tmp_path / "bin")
    query = tmp_path / "query.fasta"
    query.write_text(">q\nACGTACGTACGT\n", encoding="utf-8")
    output = tmp_path / "ptu.tsv"

    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo_root / "src")

    check = subprocess.run(
        [
            sys.executable, "-m", "plasrip.cli", "database", "--check", str(db), "--json"
        ],
        cwd=repo_root,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )
    assert '"PTUs": 1' in check.stdout
    assert '"PGGs": 1' in check.stdout

    subprocess.run(
        [
            sys.executable, "-m", "plasrip.cli", "assign",
            "-i", str(query), "--database", str(db),
            "-o", str(output), "--mash", str(fake_mash), "--quiet",
        ],
        cwd=repo_root,
        env=env,
        check=True,
        capture_output=True,
        text=True,
    )

    result = pd.read_csv(output, sep="\t", dtype=str, keep_default_na=False)
    assert result.loc[0, "PLASMID"] == "query"
    assert result.loc[0, "PTU"] == "PTU000001"
    assert result.loc[0, "PGG_ID"] == "PGG-000001"
    assert result.loc[0, "PGG_assignment_status"] == "ASSIGNED_STABLE_PGG"
    assert result.loc[0, "how_assigned"] == "all_members"

    hits = pd.read_csv(str(output) + ".hits.tsv", sep="\t", dtype=str, keep_default_na=False)
    assert hits.loc[0, "PTU_ID"] == "PTU000001"
    assert hits.loc[0, "PGG_ID"] == "PGG-000001"


def test_step9_annotation_summary_normalizes_and_merges_all_layers(tmp_path: Path):
    fasta = tmp_path / "p1.fasta"
    fasta.write_text(">p1\nACGT\n", encoding="utf-8")

    amr = pd.DataFrame([
        {"PLASMID": fasta.name, "GENE": "blaX"},
    ])
    replicon = pd.DataFrame([
        {"PLASMID": fasta.name, "GENE": "IncX"},
    ])
    mobility = pd.DataFrame([
        {"PLASMID": fasta.stem, "mobility": "pCONJ"},
    ])

    out = annotate.build_annotation_summary(
        [str(fasta)], amr_df=amr, replicon_df=replicon, mobility_df=mobility
    )
    assert list(out.columns) == ANNOTATION_SUMMARY_COLUMNS
    assert out.loc[0, "PLASMID"] == "p1"
    assert int(out.loc[0, "AMR_hits"]) == 1
    assert out.loc[0, "AMR_genes"] == "blaX"
    assert int(out.loc[0, "REPLICON_hits"]) == 1
    assert out.loc[0, "REPLICON_genes"] == "IncX"
    assert out.loc[0, "MOB"] == "pCONJ"


def test_step9_run_annotation_orchestration_with_mocked_external_tools(tmp_path: Path, monkeypatch):
    fasta = tmp_path / "p1.fasta"
    fasta.write_text(">p1\nACGTACGT\n", encoding="utf-8")

    monkeypatch.setattr(annotate, "check_dependencies", lambda **kwargs: None)
    monkeypatch.setattr(
        annotate,
        "run_amr_pipeline",
        lambda **kwargs: pd.DataFrame([{"PLASMID": fasta.name, "GENE": "blaX"}]),
    )
    monkeypatch.setattr(
        annotate,
        "run_replicon_typing",
        lambda **kwargs: pd.DataFrame([{"PLASMID": fasta.name, "GENE": "IncX"}]),
    )
    monkeypatch.setattr(
        annotate,
        "run_mobility",
        lambda **kwargs: pd.DataFrame([{"PLASMID": fasta.stem, "mobility": "pCONJ"}]),
    )

    result = annotate.run_annotation(
        input_fastas=[str(fasta)],
        out_dir=str(tmp_path / "annotation"),
        with_amr=True,
        with_replicon=True,
        with_mobility=True,
        quiet=True,
        overwrite=True,
    )
    summary = result["annotation_summary"]
    assert len(summary) == 1
    assert summary.loc[0, "PLASMID"] == "p1"
    assert int(summary.loc[0, "AMR_hits"]) == 1
    assert int(summary.loc[0, "REPLICON_hits"]) == 1
    assert summary.loc[0, "MOB"] == "pCONJ"

    written = pd.read_csv(
        tmp_path / "annotation" / "annotation_summary.tsv",
        sep="\t", dtype=str, keep_default_na=False,
    )
    assert written.loc[0, "PLASMID"] == "p1"
    assert written.loc[0, "MOB"] == "pCONJ"
