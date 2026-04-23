from __future__ import annotations

import json

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
    for cmd in ["predict", "assign", "core", "annotate", "run", "version", "models"]:
        assert cmd in help_text
