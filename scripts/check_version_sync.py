from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = (ROOT / "VERSION").read_text(encoding="utf-8").strip()

CHECKS = {
    "pyproject.toml": r'version = "([^"]+)"',
    "src/plasrip/__init__.py": r'__version__ = "([^"]+)"',
    "src/plasrip/cli.py": r'__version__ = "([^"]+)"',
    "conda.recipe/meta.yaml": r'\{% set version = "([^"]+)" %\}',
}


def extract_version(path: Path, pattern: str) -> str:
    text = path.read_text(encoding="utf-8")
    match = re.search(pattern, text)
    if not match:
        raise RuntimeError(f"Could not find version in {path}")
    return match.group(1)


def main() -> int:
    print(f"Expected version: {EXPECTED}")
    failed = False
    for rel_path, pattern in CHECKS.items():
        value = extract_version(ROOT / rel_path, pattern)
        status = "OK" if value == EXPECTED else "MISMATCH"
        print(f"[{status}] {rel_path}: {value}")
        if value != EXPECTED:
            failed = True
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
