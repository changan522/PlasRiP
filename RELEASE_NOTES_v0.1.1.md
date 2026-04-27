# PlasRiP v0.1.1

PlasRiP v0.1.1 is a maintenance and usability release focused on improving deployment behavior, host-range model discovery, resource documentation, and command-line error reporting.

## Highlights

- Automatically resolves the bundled host-range model directory when `--model-dir` is not provided.
- Supports `PLASRIP_MODEL_DIR` as an optional override for custom host-range model locations.
- Improves error messages for missing FASTA files, unexpanded wildcard paths, empty FASTA files, invalid FASTA files, and missing/incomplete host-range model directories.
- Documents the external PTU sketch library distributed as a GitHub Release asset.
- Updates README documentation with a workflow overview figure, installation instructions, quick-start examples, output descriptions, and troubleshooting notes.

## Recommended installation

```bash
conda create -n plasrip python=3.11 -y
conda activate plasrip
conda install -c changan522 -c conda-forge -c bioconda plasrip=0.1.1
```

Verify the installation:

```bash
plasrip --help
plasrip version
plasrip models --json
```

## PTU sketch library

PTU assignment requires an external PTU sketch library. Download the release asset:

```bash
mkdir -p ~/plasrip_resources
cd ~/plasrip_resources

wget -O PlasRiP_PTU_sketch_library_v0.1.1.tar.gz \
  https://github.com/changan522/PlasRiP/releases/download/v0.1.1/PlasRiP_PTU_sketch_library_v0.1.1.tar.gz

wget -O PlasRiP_PTU_sketch_library_v0.1.1.sha256 \
  https://github.com/changan522/PlasRiP/releases/download/v0.1.1/PlasRiP_PTU_sketch_library_v0.1.1.sha256

sha256sum -c PlasRiP_PTU_sketch_library_v0.1.1.sha256
tar -xzf PlasRiP_PTU_sketch_library_v0.1.1.tar.gz

PTU_REP_DIR=$(tar -tzf PlasRiP_PTU_sketch_library_v0.1.1.tar.gz | head -1 | cut -d/ -f1)
```

## Main usage example

```bash
plasrip run \
  -i examples/plasmid_genomes/pCZ2103-1.fasta \
  --rep-dir ~/plasrip_resources/${PTU_REP_DIR} \
  -o run_out \
  -t 8
```

## Notes

- The packaged host-range model is resolved automatically from `plasrip/data/hostrange`.
- PTU assignment requires the external PTU sketch library through `--rep-dir`.
- External tools such as Mash, ABRicate, Prodigal, and MacSyFinder are expected to be installed through the recommended conda workflow.
