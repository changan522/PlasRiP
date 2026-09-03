# PlasRiP v0.2.0

PlasRiP v0.2.0 introduces a versioned external reference-database architecture for PTU assignment and adds higher-order PGG contextualization while preserving the established host-range prediction and annotation workflows.

## Highlights

- Added support for the external **PlasRiP Reference Database v2.0**:
  <https://github.com/changan522/PlasRiP/releases/tag/database-v2.0>
- Replaced the legacy PTU sketch-directory interface with `--database` and `PLASRIP_DATABASE`.
- Added centralized database loading and validation through `plasrip database`.
- Added formal PTU identifiers (`PTU000001` style) and PGG identifiers (`PGG-000001` style).
- PTU assignment now reports `PGG_ID` and `PGG_assignment_status` when available.
- Query Mash sketches are validated against database-manifest parameters.
- Separated model-predicted host spectra (`Predicted_*`) from observed PTU host metadata (`PTU_observed_*`).
- Removed the uncalibrated `Predicted_host_confidence` heuristic after server validation; raw model scores remain available.
- Simplified the user-facing `summary.tsv` to 25 high-value columns while retaining detailed diagnostics in dedicated output files.
- Removed legacy bundled PTU metadata (`PTU-metadata.csv`, `PTUmember.csv`, and `PTU_info.tsv`) from the software package.

## Reference Database v2.0

The companion database contains:

- 64,391 reference plasmids
- 52,591 plasmids assigned to 5,975 PTUs
- 11,800 reference plasmids without stable PTU assignment
- 4,152 PTUs organized into 185 PGGs

Reference Mash parameters:

- Mash version 2.3
- k-mer size 21
- sketch size 1000
- seed 42

Download and documentation:

<https://github.com/changan522/PlasRiP/releases/tag/database-v2.0>

## Compatibility note

v0.2.0 changes the PTU reference-resource interface. Workflows written for v0.1.x that use `--rep-dir` or package-internal PTU metadata should migrate to the external v2.0 database and `--database` interface.

The software and database are versioned independently:

- software tag: `v0.2.0`
- database tag: `database-v2.0`
