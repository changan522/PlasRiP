# PlasRiP v0.3.0 conda release guide

## 1. Verify version synchronization

```bash
python scripts/check_version_sync.py
```

The canonical version must be `0.3.0` and should match `VERSION`, `pyproject.toml`, `src/plasrip/__init__.py`, and `conda.recipe/meta.yaml`.

## 2. Run tests

```bash
pytest -q
```

## 3. Build the conda package

```bash
conda install -n base -c conda-forge conda-build anaconda-client
conda build conda.recipe
```

Print the generated artifact path:

```bash
conda build conda.recipe --output
```

## 4. Test a local conda installation

```bash
conda create -n plasrip-test python=3.10 -y
conda activate plasrip-test
conda install --use-local -c conda-forge -c bioconda plasrip=0.3.0
plasrip version
plasrip models --json
```

The external PlasRiP Reference Database v3.0 is distributed separately and is not bundled in the conda package.

## 5. Upload to Anaconda.org

```bash
anaconda login
anaconda upload $(conda build conda.recipe --output)
```

## 6. GitHub release convention

Software tag: `v0.3.0`  
Software release title: `PlasRiP v0.3.0`  
Database tag: `database-v3.0`

The software and database are versioned independently.
