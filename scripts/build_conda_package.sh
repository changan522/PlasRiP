#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

python scripts/check_version_sync.py
conda install -n base -y -c conda-forge conda-build anaconda-client
conda build conda.recipe

echo
echo "Build finished for PlasRip v$(cat VERSION)."
echo "To locate the package, run:"
echo "  conda build $ROOT_DIR/conda.recipe --output"
