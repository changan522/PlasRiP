#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

python scripts/check_version_sync.py
PKG_PATH="$(conda build conda.recipe --output)"
echo "Uploading package: $PKG_PATH"
anaconda upload "$PKG_PATH"
