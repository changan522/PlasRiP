# PlasRiP v0.2.1

PlasRiP v0.2.1 is a packaging and release-consistency update to v0.2.0. It adds a fully supported Conda installation path and explicitly pins Mash 2.3 to match the companion PlasRiP Reference Database v2.0.

## Highlights

- PlasRiP v0.2.1 can be installed directly from the `changan522` Conda channel.
- The Conda/runtime Mash dependency is explicitly pinned to **Mash 2.3**.
- GitHub, Conda, and database installation documentation are synchronized.
- The companion reference database remains **PlasRiP Reference Database v2.0**.
- No scientific algorithms or public output schemas were changed from v0.2.0.

## Recommended installation

```bash
conda create -n plasrip \
  -c changan522 \
  -c conda-forge \
  -c bioconda \
  plasrip=0.2.1 -y

conda activate plasrip
plasrip --version
```

The Conda installation resolves the required runtime tools, including Mash 2.3, ABRicate, Prodigal, and MacSyFinder.

## Reference database

PlasRiP v0.2.1 uses **PlasRiP Reference Database v2.0**:

<https://github.com/changan522/PlasRiP/releases/tag/database-v2.0>

The v2.0 database contains:

- 64,391 reference plasmids
- 52,591 plasmids assigned to 5,975 plasmid taxonomic units (PTUs)
- 11,800 reference plasmids without stable PTU assignment
- 4,152 PTUs organized into 185 plasmid genome groups (PGGs)

Reference Mash parameters:

- Mash version 2.3
- k-mer size 21
- sketch size 1000
- seed 42

## Alternative wheel installation

The wheel attached to this release can be installed with:

```bash
pip install plasrip-0.2.1-py3-none-any.whl
```

When installing from the wheel, external bioinformatics tools must be installed separately.

## Compatibility

- Software: **PlasRiP v0.2.1**
- Reference database: **v2.0**
- Database schema: **2.0**
- Mash: **2.3**

No host-range prediction, PTU assignment, PGG contextualization, AMR, replicon, mobility, or summary-schema algorithms were changed from v0.2.0.

## Full changelog

https://github.com/changan522/PlasRiP/compare/v0.2.0...v0.2.1

## Citation

Citation information will be updated upon publication of the associated manuscript.
