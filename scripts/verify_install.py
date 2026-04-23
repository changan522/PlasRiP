#!/usr/bin/env python3
"""Minimal post-install verification for PlasRip."""

from __future__ import annotations

import contextlib
import json
import sys
from io import StringIO
from pathlib import Path

import plasrip
from plasrip import cli


def main() -> int:
    print(f"[OK] Imported plasrip {plasrip.__version__}")

    rc = cli.main(["version"])
    if rc != 0:
        print("[FAIL] `plasrip version` failed", file=sys.stderr)
        return rc
    print("[OK] `plasrip version` completed")

    buf = StringIO()
    with contextlib.redirect_stdout(buf):
        rc = cli.main(["models", "--json"])
    if rc != 0:
        print("[FAIL] `plasrip models --json` failed", file=sys.stderr)
        return rc

    payload = json.loads(buf.getvalue())
    print("[OK] `plasrip models --json` completed")
    print("[INFO] Dependency report keys:", ", ".join(sorted(payload.keys())))

    repo_root = Path(__file__).resolve().parents[1]
    plasmid_dir = repo_root / "examples" / "plasmid_genomes"
    if plasmid_dir.is_dir():
        print(f"[OK] Example plasmid genome directory exists: {plasmid_dir}")
    else:
        print(f"[WARN] Example plasmid genome directory is missing: {plasmid_dir}")

    print("\nMinimal installation verification finished.")
    print("Next step: place plasmid FASTA files under examples/plasmid_genomes/ and run a real command.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
