# PlasRip conda release guide

## 1. Verify version synchronization

```bash
python scripts/check_version_sync.py
```

The canonical version is stored in `VERSION`. For this release it must be `0.1.0`, and it should match:

- `pyproject.toml`
- `src/plasrip/__init__.py`
- `src/plasrip/cli.py`
- `conda.recipe/meta.yaml`

## 2. Create the build environment

```bash
conda env create -f environment.yml
conda activate plasrip
```

## 3. Run pre-release checks

```bash
bash scripts/release_prep.sh
```

## 4. Build the conda package locally

```bash
bash scripts/build_conda_package.sh
```

To print the exact output artifact path:

```bash
conda build conda.recipe --output
```

## 5. Test installation from the locally built package

```bash
bash scripts/install_local_conda_package.sh plasrip-test
conda activate plasrip-test
plasrip --help
plasrip version
python scripts/verify_install.py
```

## 6. Upload to Anaconda.org

```bash
anaconda login
bash scripts/upload_anaconda.sh
```

## 7. Tag the release on GitHub

```bash
git add .
git commit -m "release: v0.1.0"
git tag -a v0.1.0 -m "PlasRip v0.1.0"
git push origin main
git push origin v0.1.0
```

## 8. Recommended release checklist

- Confirm `python scripts/check_version_sync.py` passes
- Confirm all packaged data assets are present under `src/plasrip/data/`
- Run smoke tests defined in `conda.recipe/meta.yaml`
- Remove local build outputs before creating a release archive
- Use an annotated Git tag with the `v` prefix
- Paste `RELEASE_NOTES_v0.1.0.md` into the GitHub release description

## Notes

- This recipe is intended for Linux or WSL environments.
- Windows builds are intentionally disabled in `conda.recipe/meta.yaml`.
- The package build uses `pip install . --no-deps` because dependencies are managed by conda.
