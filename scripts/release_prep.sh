#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

python scripts/check_version_sync.py
python -m compileall src tests scripts
pytest -q tests/test_import_and_version.py tests/test_cli.py tests/test_install_smoke.py tests/test_validate.py tests/test_schemas.py

echo
echo "Release prep checks passed for version $(cat VERSION)"
