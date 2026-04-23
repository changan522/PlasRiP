PlasRip conda recipe notes
==========================

This directory contains a Linux/WSL-oriented conda recipe for building the
PlasRip package from the local source tree.

Key design choices:
- The package is declared as `noarch: python`, but Windows builds are disabled.
- Dependencies are managed by conda; the build step uses `pip install . --no-deps`.
- Smoke tests are executed during `conda build` via `meta.yaml`.
- Model/data assets under `src/plasrip/data/` are packaged via `pyproject.toml`.
