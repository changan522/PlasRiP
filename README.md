# PlasRip

**Version:** 0.1.0

PlasRip is a plasmid analysis platform centered on two core capabilities:

- host-range prediction
- PTU assignment

It also integrates one-stop plasmid annotation using established tools for:

- AMR detection
- replicon typing
- mobility prediction

## Release layout

This release tree is prepared for two publication targets:

- **GitHub** for source-code hosting, tagging, and release notes
- **Anaconda.org** for distributing the conda package `plasrip`

Key release files:

- `VERSION` — canonical release version
- `CHANGELOG.md` — project change history
- `RELEASE_NOTES_v0.1.0.md` — GitHub release body for this version
- `TAGGING_AND_RELEASE.md` — tag naming and release convention
- `ANACONDA_UPLOAD.md` — Anaconda.org upload commands
- `scripts/check_version_sync.py` — validates that all version declarations match

## Conda installation

### Option A: create the runtime environment and install from source

```bash
conda env create -f environment.yml
conda activate plasrip
pip install . --no-deps
```

For development mode:

```bash
pip install -e . --no-deps
```

### Option B: build a conda package

```bash
conda install -n base -c conda-forge conda-build anaconda-client
conda build conda.recipe
conda install --use-local plasrip
```

A release-oriented build/test guide is provided in `CONDA_RELEASE_GUIDE.md`.

## GitHub release workflow

Recommended release flow for version `0.1.0`:

```bash
python scripts/check_version_sync.py
bash scripts/release_prep.sh
git add .
git commit -m "release: v0.1.0"
git tag -a v0.1.0 -m "PlasRip v0.1.0"
git push origin main
git push origin v0.1.0
```

Then create a GitHub release titled **PlasRip v0.1.0** and paste the contents of `RELEASE_NOTES_v0.1.0.md` into the release description.

## Anaconda.org upload workflow

```bash
conda install -n base -y -c conda-forge conda-build anaconda-client
conda build conda.recipe
anaconda login
anaconda upload $(conda build conda.recipe --output)
```

For a single-command upload helper:

```bash
bash scripts/upload_anaconda.sh
```

## External tools

Some modules require external command-line tools to be available in `PATH`:

- `mash` for PTU assignment
- `abricate` for AMR and replicon typing
- `prodigal` for ORF prediction
- `macsyfinder` for mobility prediction

These are already included in the provided `environment.yml` and `conda.recipe/meta.yaml`.

## Data and model files

PTU assignment expects:

- packaged metadata files under `src/plasrip/data/`
- a user-supplied PTU `.msh` sketch directory via `--rep-dir`

Host-range prediction expects a model directory provided via `--model-dir`. The model directory can either be a structured directory or a flat directory containing files such as:

- `two_tower_best_3.2.pth`
- `host_profiles.joblib`
- `host_vocab.joblib`
- `plasmid_codon_scaler.joblib`

## Command line entry

After installation, the canonical entry point is:

```bash
plasrip
```

Available subcommands:

```bash
plasrip predict
plasrip assign
plasrip core
plasrip annotate
plasrip run
plasrip version
plasrip models
```

## Quick examples

### Host-range prediction

```bash
plasrip predict \
  -i plasmid1.fasta plasmid2.fasta \
  --model-dir /path/to/host_model \
  -o predict_out
```

### PTU assignment

```bash
plasrip assign \
  -i plasmid1.fasta plasmid2.fasta \
  -r /path/to/ptu_msh_dir \
  -o ptu.tsv
```

### Dual-core analysis

```bash
plasrip core \
  -i plasmid1.fasta plasmid2.fasta \
  --model-dir /path/to/host_model \
  --rep-dir /path/to/ptu_msh_dir \
  -o core_out
```

### Full workflow

```bash
plasrip run \
  -i plasmid1.fasta plasmid2.fasta \
  --model-dir /path/to/host_model \
  --rep-dir /path/to/ptu_msh_dir \
  -o out_dir
```

## Output overview

Typical outputs include:

- `hostrange.tsv`
- `hostrange.top_hits.tsv`
- `hostrange.taxonomy.tsv`
- `ptu.tsv`
- `ptu.tsv.hits.tsv`
- `core_summary.tsv`
- `amr.tsv`
- `plasmidfinder.tsv`
- `mobility_summary.tsv`
- `annotation_summary.tsv`
- `summary.tsv`
- `PlasRip.log`

## Legacy entry points

For backward compatibility, the package also exposes several direct script entry points such as:

- `PlasRip`
- `PlasRip-assign`
- `PlasRip-amr`
- `PlasRip-replicon`
- `PlasRip-predict`
- `PlasRip-core`
- `PlasRip-annotate`
- `PlasRip-run`

The recommended entry point is `plasrip`.
