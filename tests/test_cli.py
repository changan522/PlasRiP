from __future__ import annotations

import json

import pytest

from plasrip import cli


def test_cli_version(capsys):
    rc = cli.main(["version"])
    captured = capsys.readouterr()
    assert rc == 0
    assert "PlasRiP" in (captured.out + captured.err)


def test_cli_models_json(capsys):
    rc = cli.main(["models", "--json"])
    captured = capsys.readouterr()
    assert rc == 0
    payload = json.loads(captured.out)
    assert "plasrip_version" in payload
    assert "dependencies" in payload


def test_cli_parser_has_expected_commands():
    parser = cli.build_parser()
    help_text = parser.format_help()
    for cmd in ["predict", "assign", "core", "annotate", "run", "version", "models", "database"]:
        assert cmd in help_text


def test_assign_uses_database_and_consistent_ptu_threshold_name():
    parser = cli.build_parser()
    args = parser.parse_args([
        "assign",
        "-i", "query.fasta",
        "--database", "/db",
        "-o", "ptu.tsv",
        "--ptu-threshold", "0.04",
    ])
    assert args.database == "/db"
    assert args.ptu_threshold == 0.04
    assert not hasattr(args, "rep_dir")


def test_legacy_rep_dir_is_not_accepted():
    parser = cli.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args([
            "assign",
            "-i", "query.fasta",
            "--rep-dir", "/old/sketches",
            "-o", "ptu.tsv",
        ])


def test_database_check_explicit_path_json(monkeypatch, capsys):
    seen = {}

    def fake_report(database, deep=False):
        seen["database"] = database
        seen["deep"] = deep
        return {"root": database, "database_version": "v2.0", "PTUs": 5975, "PGGs": 185}

    monkeypatch.setattr(cli, "_database_report", fake_report)
    rc = cli.main(["database", "--check", "/db", "--deep", "--json"])
    captured = capsys.readouterr()
    assert rc == 0
    assert seen == {"database": "/db", "deep": True}
    payload = json.loads(captured.out)
    assert payload["database_version"] == "v2.0"


def test_database_check_can_use_environment_resolution(monkeypatch, capsys):
    seen = {}

    def fake_report(database, deep=False):
        seen["database"] = database
        seen["deep"] = deep
        return {"root": "/env/db", "database_version": "v2.0", "PTUs": 5975, "PGGs": 185}

    monkeypatch.setattr(cli, "_database_report", fake_report)
    rc = cli.main(["database", "--check", "--json"])
    captured = capsys.readouterr()
    assert rc == 0
    assert seen == {"database": None, "deep": False}
    payload = json.loads(captured.out)
    assert payload["root"] == "/env/db"


def test_database_command_requires_check_flag():
    parser = cli.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["database"])
