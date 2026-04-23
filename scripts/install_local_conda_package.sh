#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

ENV_NAME="${1:-plasrip-test}"
python scripts/check_version_sync.py
PKG_PATH="$(conda build "$ROOT_DIR/conda.recipe" --output)"

conda create -n "$ENV_NAME" -y -c conda-forge -c bioconda python=3.10
conda install -n "$ENV_NAME" -y -c conda-forge -c bioconda "$PKG_PATH"

echo
echo "Installed package into environment: $ENV_NAME"
echo "Activate with: conda activate $ENV_NAME"
