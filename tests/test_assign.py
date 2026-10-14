from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from plasrip import assign
from plasrip.exceptions import PlasRipConfigError


def _hit(mpc: str, dist: float, rep: str | None = None):
    return {
        "plasmid": "query.fasta",
        "MPC_ID": mpc,
        "Representative": rep or f"rep_{mpc}",
        "dist": dist,
    }


def test_reference_batches_preserve_order_and_split_by_size():
    refs = ["a" * 6, "b" * 6, "c" * 6]
    batches = assign._reference_batches(refs, max_command_chars=14)
    assert batches == [[refs[0], refs[1]], [refs[2]]]
    assert [x for batch in batches for x in batch] == refs


def test_mash_dist_multi_combines_batch_outputs(monkeypatch):
    calls = []

    def fake_batches(ref_paths):
        assert ref_paths == ["r1.msh", "r2.msh", "r3.msh"]
        return [["r1.msh", "r2.msh"], ["r3.msh"]]

    def fake_run(cmd, check, capture_output, text):
        calls.append(cmd)
        if "r3.msh" in cmd:
            stdout = "query\tref3\t0.03\t0\t1/1\n"
        else:
            stdout = "query\tref1\t0.01\t0\t1/1\nquery\tref2\t0.02\t0\t1/1\n"
        return subprocess.CompletedProcess(cmd, 0, stdout=stdout, stderr="")

    monkeypatch.setattr(assign, "_reference_batches", fake_batches)
    monkeypatch.setattr(assign.subprocess, "run", fake_run)

    lines = assign.mash_dist_multi(
        "query.msh", ["r1.msh", "r2.msh", "r3.msh"], "mash", 2
    )
    assert len(calls) == 2
    assert lines == [
        "query\tref1\t0.01\t0\t1/1",
        "query\tref2\t0.02\t0\t1/1",
        "query\tref3\t0.03\t0\t1/1",
    ]


def test_assign_mpc_no_hit_and_failed():
    no_hit = assign.assign_mpc(
        "query.fasta", "query.msh", [], {}, 0.05, "mash", 1
    )
    assert no_hit["mpc_assigned"] == "NA"
    assert no_hit["how_assigned"] == "no_hit"

    failed = assign.assign_mpc(
        "query.fasta", None, [], {}, 0.05, "mash", 1
    )
    assert failed["mpc_assigned"] == "NA"
    assert failed["how_assigned"] == "failed"


def test_assign_mpc_full_member_path_is_preserved(monkeypatch):
    monkeypatch.setattr(
        assign,
        "full_member_check",
        lambda *args, **kwargs: (True, 0.02, 0.01, 2, 1.0, 2, 0.01),
    )
    result = assign.assign_mpc(
        "query.fasta",
        "query.msh",
        [_hit("MPC000001", 0.01)],
        {"MPC000001": ["a.msh", "b.msh"]},
        0.05,
        "mash",
        1,
    )
    assert result["mpc_assigned"] == "MPC000001"
    assert result["how_assigned"] == "all_members"
    assert result["all_member_n"] == 2


def test_assign_mpc_single_candidate_rep_fallback_is_preserved(monkeypatch):
    monkeypatch.setattr(
        assign,
        "full_member_check",
        lambda *args, **kwargs: (False, 0.08, 0.05, 2, 0.5, 1, 0.02),
    )
    result = assign.assign_mpc(
        "query.fasta",
        "query.msh",
        [_hit("MPC000001", 0.01)],
        {"MPC000001": ["a.msh", "b.msh"]},
        0.05,
        "mash",
        1,
    )
    assert result["mpc_assigned"] == "MPC000001"
    assert result["how_assigned"] == "single_mpc_by_rep"


def test_assign_mpc_network_support_path_is_preserved(monkeypatch):
    stats = {
        "MPC000001": (False, 0.08, 0.04, 5, 0.8, 4, 0.025),
        "MPC000002": (False, 0.09, 0.05, 5, 0.4, 2, 0.03),
    }

    def fake_full(query_msh, member_paths, threshold, mash_cmd, threads):
        return stats[member_paths[0]]

    monkeypatch.setattr(assign, "full_member_check", fake_full)
    result = assign.assign_mpc(
        "query.fasta",
        "query.msh",
        [_hit("MPC000001", 0.02), _hit("MPC000002", 0.025)],
        {
            "MPC000001": ["MPC000001"],
            "MPC000002": ["MPC000002"],
        },
        0.05,
        "mash",
        1,
    )
    assert result["mpc_assigned"] == "MPC000001"
    assert result["how_assigned"] == "network_support"
    assert result["network_coverage"] == 0.8


def test_assign_mpc_closest_rep_and_ambiguous_paths_are_preserved(monkeypatch):
    monkeypatch.setattr(
        assign,
        "full_member_check",
        lambda *args, **kwargs: (False, 0.09, 0.08, 2, 0.0, 0, None),
    )
    members = {
        "MPC000001": ["a.msh", "b.msh"],
        "MPC000002": ["c.msh", "d.msh"],
    }

    closest = assign.assign_mpc(
        "query.fasta",
        "query.msh",
        [_hit("MPC000001", 0.01), _hit("MPC000002", 0.03)],
        members,
        0.05,
        "mash",
        1,
    )
    assert closest["mpc_assigned"] == "MPC000001"
    assert closest["how_assigned"] == "closest_rep"

    ambiguous = assign.assign_mpc(
        "query.fasta",
        "query.msh",
        [_hit("MPC000001", 0.01), _hit("MPC000002", 0.015)],
        members,
        0.05,
        "mash",
        1,
    )
    assert ambiguous["mpc_assigned"] == "NA"
    assert ambiguous["how_assigned"] == "ambiguous"


def test_direct_assign_cli_uses_step2_mpc_threshold_name():
    args = assign.parse_args([
        "-i", "query.fasta",
        "--database", "/db",
        "-o", "mpc.tsv",
        "--mpc-threshold", "0.04",
    ])
    assert args.mpc_threshold == 0.04
    with pytest.raises(SystemExit):
        assign.parse_args([
            "-i", "query.fasta",
            "--database", "/db",
            "-o", "mpc.tsv",
            "--threshold", "0.04",
        ])


def test_run_assignment_rejects_invalid_cutoff_before_external_mash(monkeypatch):
    class FakeDB:
        mpc_threshold = 0.05
        mash_k = 21
        mash_sketch_size = 1000
        mash_seed = 42
        root = "/db"
        database_version = "v3.0"
        mash_index_dir = "/db/mash_index"

    monkeypatch.setattr(assign, "load_database", lambda *a, **k: FakeDB())
    monkeypatch.setattr(assign, "validate_query_mash_parameters", lambda *a, **k: (21, 1000))
    monkeypatch.setattr(assign, "validate_mash_runtime", lambda *a, **k: "2.3")
    monkeypatch.setattr(assign, "validate_reference_sketch_sample", lambda *a, **k: [])

    with pytest.raises(PlasRipConfigError, match="within \\[0, 1\\]"):
        assign.run_mpc_assignment(
            ["does-not-matter.fasta"], "/db", "out.tsv", threshold=1.1, show_progress=False
        )


def test_run_mpc_assignment_uses_database_indexes_and_writes_summary(tmp_path, monkeypatch):
    fasta = tmp_path / "query.fasta"
    fasta.write_text(">query\nACGTACGTACGT\n", encoding="utf-8")
    output = tmp_path / "mpc.tsv"

    class FakeDB:
        mpc_threshold = 0.05
        mash_k = 21
        mash_sketch_size = 1000
        mash_seed = 42
        root = tmp_path / "db"
        database_version = "v3.0"
        mash_index_dir = tmp_path / "db" / "mash_index"

    fake_db = FakeDB()
    monkeypatch.setattr(assign, "validate_query_mash_parameters", lambda *a, **k: (21, 1000))
    monkeypatch.setattr(assign, "validate_mash_runtime", lambda *a, **k: "2.3")
    monkeypatch.setattr(assign, "validate_reference_sketch_sample", lambda *a, **k: [])
    monkeypatch.setattr(assign, "ensure_sketch", lambda *a, **k: str(tmp_path / "query.msh"))
    monkeypatch.setattr(
        assign,
        "mash_dist_multi",
        lambda *a, **k: [
            "query.fasta\tplasmid_fasta/ref1.fasta\t0.01\t0\t1000/1000"
        ],
    )
    monkeypatch.setattr(
        assign,
        "full_member_check",
        lambda *a, **k: (True, 0.02, 0.015, 2, 1.0, 2, 0.015),
    )

    df = assign.run_mpc_assignment(
        [str(fasta)],
        str(fake_db.root),
        str(output),
        show_progress=False,
        _database=fake_db,
        _assignment_indexes=(
            {"MPC000001": ["ref1.msh", "ref2.msh"]},
            ["ref1.msh"],
            {"plasmid_fasta/ref1.fasta": ("MPC000001", "ref1")},
        ),
        _mpc_meta_lookup={
            "MPC000001": {
                "MPC_ID": "MPC000001",
                "PGG_ID": "PGG-000001",
                "PGG_assignment_status": "ASSIGNED_STABLE_PGG",
            }
        },
    )

    assert output.is_file()
    assert list(df.columns) == assign.MPC_COLUMNS
    assert df.loc[0, "MPC"] == "MPC000001"
    assert df.loc[0, "PGG_ID"] == "PGG-000001"
    assert df.loc[0, "PGG_assignment_status"] == "ASSIGNED_STABLE_PGG"
    assert df.loc[0, "how_assigned"] == "all_members"

    hits_path = tmp_path / "mpc.tsv.hits.tsv"
    assert hits_path.is_file()
    hits = __import__("pandas").read_csv(hits_path, sep="\t", keep_default_na=False)
    assert list(hits.columns) == assign.MPC_HITS_COLUMNS
    assert hits.loc[0, "MPC_ID"] == "MPC000001"
    assert hits.loc[0, "PGG_ID"] == "PGG-000001"
    assert hits.loc[0, "PGG_assignment_status"] == "ASSIGNED_STABLE_PGG"


def test_ensure_sketch_builds_with_explicit_database_seed_k_and_size(tmp_path, monkeypatch):
    fasta = tmp_path / "query.fasta"
    fasta.write_text(">q\nACGTACGT\n", encoding="utf-8")
    commands = []

    class FakeDB:
        database_version = "v3.0"

    def fake_run(cmd, check, capture_output, text):
        commands.append(cmd)
        # emulate Mash creating the requested output sketch
        prefix = cmd[cmd.index("-o") + 1]
        Path(prefix + ".msh").write_bytes(b"fake")
        return subprocess.CompletedProcess(cmd, 0, stdout="", stderr="")

    class FakeInfo:
        pass

    monkeypatch.setattr(assign.subprocess, "run", fake_run)
    monkeypatch.setattr(assign, "inspect_mash_sketch", lambda *a, **k: FakeInfo())
    monkeypatch.setattr(assign, "validate_mash_sketch_compatibility", lambda *a, **k: None)

    msh = assign.ensure_sketch(
        str(fasta), "mash", 4, mash_k=21, mash_s=1000, mash_seed=42, database=FakeDB()
    )
    assert msh.endswith("query.msh")
    cmd = commands[0]
    assert cmd[:2] == ["mash", "sketch"]
    assert ["-S", "42"] == cmd[cmd.index("-S"):cmd.index("-S")+2]
    assert ["-k", "21"] == cmd[cmd.index("-k"):cmd.index("-k")+2]
    assert ["-s", "1000"] == cmd[cmd.index("-s"):cmd.index("-s")+2]


def test_ensure_sketch_rejects_incompatible_existing_query_sketch(tmp_path, monkeypatch):
    fasta = tmp_path / "query.fasta"
    fasta.write_text(">q\nACGT\n", encoding="utf-8")
    msh = tmp_path / "query.msh"
    msh.write_bytes(b"old")

    monkeypatch.setattr(assign, "inspect_mash_sketch", lambda *a, **k: object())

    def reject(*args, **kwargs):
        raise PlasRipConfigError("Existing query sketch is incompatible")

    monkeypatch.setattr(assign, "validate_mash_sketch_compatibility", reject)
    with pytest.raises(PlasRipConfigError, match="incompatible"):
        assign.ensure_sketch(
            str(fasta), "mash", 1, mash_k=21, mash_s=1000, mash_seed=42, database=object()
        )


def test_step5_public_output_preserves_unassigned_pgg_status_for_assigned_mpc():
    raw = [{
        "query": "query.fasta",
        "n_mpc_hits": 1,
        "mpc_assigned": "MPC000003",
        "PGG_ID": "",
        "PGG_assignment_status": "UNASSIGNED_OUTSIDE_CONSENSUS_SPACE",
        "how_assigned": "closest_rep",
        "closest_rep_dist": 0.01,
        "closest_mpc": "MPC000003",
        "closest_rep": "ref3",
        "second_mpc": "NA",
        "second_rep_dist": "NA",
        "rep_dist_gap": "NA",
        "all_member_mpc": "NA",
        "all_member_max_dist": "NA",
        "all_member_avg_dist": "NA",
        "all_member_n": "NA",
        "network_mpc": "NA",
        "network_coverage": "NA",
        "network_mean_dist_leq_threshold": "NA",
        "network_n_within_threshold": "NA",
    }]
    out = assign._format_assignment_output(raw)
    assert out.loc[0, "MPC"] == "MPC000003"
    assert out.loc[0, "PGG_ID"] == ""
    assert out.loc[0, "PGG_assignment_status"] == "UNASSIGNED_OUTSIDE_CONSENSUS_SPACE"


def test_step5_public_output_marks_pgg_not_applicable_when_no_mpc():
    raw = [{
        "query": "query.fasta",
        "n_mpc_hits": 0,
        "mpc_assigned": "NA",
        "how_assigned": "no_hit",
        "closest_rep_dist": "NA",
        "closest_mpc": "NA",
        "closest_rep": "NA",
        "second_mpc": "NA",
        "second_rep_dist": "NA",
        "rep_dist_gap": "NA",
        "all_member_mpc": "NA",
        "all_member_max_dist": "NA",
        "all_member_avg_dist": "NA",
        "all_member_n": "NA",
        "network_mpc": "NA",
        "network_coverage": "NA",
        "network_mean_dist_leq_threshold": "NA",
        "network_n_within_threshold": "NA",
    }]
    out = assign._format_assignment_output(raw)
    assert out.loc[0, "MPC"] == "NA"
    assert out.loc[0, "PGG_ID"] == "NA"
    assert out.loc[0, "PGG_assignment_status"] == "NOT_APPLICABLE_NO_MPC"


def test_step5_hits_output_contains_pgg_status_and_frozen_column_order():
    hits = [{
        "plasmid": "query.fasta",
        "MPC_ID": "MPC000002",
        "PGG_ID": "",
        "PGG_assignment_status": "UNASSIGNED_NOT_ELIGIBLE",
        "Representative": "ref2",
        "dist": 0.02,
    }]
    out = assign._format_hits_output(hits)
    assert list(out.columns) == assign.MPC_HITS_COLUMNS
    assert out.loc[0, "PLASMID"] == "query"
    assert out.loc[0, "PGG_assignment_status"] == "UNASSIGNED_NOT_ELIGIBLE"


def test_step5_empty_hits_output_still_has_public_schema():
    out = assign._format_hits_output([])
    assert out.empty
    assert list(out.columns) == assign.MPC_HITS_COLUMNS
