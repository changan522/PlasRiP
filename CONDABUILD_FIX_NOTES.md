# Conda packaging fixes for PlasRiP 0.1.0

This package only changes conda/release metadata and does not alter runtime Python source code under `src/plasrip`.

## Changes made

1. Removed the duplicate `build.script` entry from `conda.recipe/meta.yaml`.
   - The recipe now uses `conda.recipe/build.sh` as the single build entrypoint.
2. Restricted conda build/install targets to Linux.
   - `skip: true  # [win or osx]`
3. Pinned Python for conda packaging to a tested range.
   - `>=3.10,<3.12`
   - Added `conda.recipe/conda_build_config.yaml` with Python 3.11 and NumPy 1.26.
4. Updated repository URLs from the old placeholder repo to:
   - `https://github.com/changan522/PlasRiP`

## Recommended commands

```bash
conda install -n base -c conda-forge conda-build anaconda-client
conda build conda.recipe
conda install --use-local plasrip
plasrip --help
plasrip version
plasrip models --json
```
