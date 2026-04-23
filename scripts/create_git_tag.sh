#!/bin/bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

VERSION="${1:-$(cat VERSION)}"
TAG="v${VERSION#v}"

python scripts/check_version_sync.py
git tag -a "$TAG" -m "PlasRip $TAG"
echo "Created annotated tag: $TAG"
