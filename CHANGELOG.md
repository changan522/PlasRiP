# Changelog

## v0.3.0 — 2026-10-14

- Rename the Mash-based reference clusters from PTU to MPC throughout the current software interface.
- Read MPC fields, metadata filenames and assignment statuses from Reference Database v3.0, schema 2.0.
- Rename assignment outputs to `mpc.tsv` and `mpc.tsv.hits.tsv`, with MPC fields in combined reports.
- Rename CLI options and Python API parameters to MPC terminology. This is an interface migration from v0.2.1.
- Preserve assignment algorithms, Mash reference paths, model weights and application score-gap default (0.20).
- Update wheel/build metadata and conda recipe to v0.3.0. Historical release notes retain their original terminology.


All notable changes to PlasRiP are documented in this file.


## [0.2.1] - 2026-09-03

### Changed
- Added first-class Conda installation instructions for the `changan522` channel.
- Pinned the Conda/runtime Mash dependency to version 2.3 to match PlasRiP Reference Database v2.0.
- Updated release and installation documentation so the GitHub, Conda, and database distributions describe the same supported configuration.

### Compatibility
- No host-range prediction, PTU assignment, PGG contextualization, annotation, or output-schema algorithms were changed from v0.2.0.
- Continues to use PlasRiP Reference Database v2.0 (schema 2.0).

## [0.2.0] - 2026-09-03

### Added
- External versioned PlasRiP reference-database architecture (schema 2.0).
- `plasrip database` command for reference-database validation.
- `PLASRIP_DATABASE` environment-variable support.
- PGG contextualization for assigned PTUs (`PGG_ID` and `PGG_assignment_status`).
- Database-backed PTU/PGG metadata integration in core outputs.
- Mash runtime compatibility checks against database-manifest parameters.
- Public Python reference-database API (`PlasRipReferenceDatabase`, `load_reference_database`, and `validate_reference_database`).

### Changed
- Replaced `--rep-dir` with `--database` for PTU-related workflows.
- PTU assignment now reads `plasmid_reference_index.tsv`, `PTU_metadata.tsv`, and `PGG_metadata.tsv` from the external database.
- Formal PTU identifiers use `PTU000001` style; PGG identifiers use `PGG-000001` style.
- Core outputs explicitly distinguish model-predicted host spectra (`Predicted_*`) from observed PTU host context (`PTU_observed_*`).
- Simplified `summary.tsv` from 54 to 25 user-facing columns; detailed diagnostics remain in `hostrange.tsv`, `ptu.tsv`, and `core_summary.tsv`.
- Query Mash parameters are read/validated against `database_manifest.tsv`.

### Removed
- Uncalibrated `Predicted_host_confidence` heuristic from public outputs.
- Bundled legacy PTU files `PTU-metadata.csv`, `PTUmember.csv`, and `PTU_info.tsv`.
- Legacy `--rep-dir` PTU-resource interface.

### Validation
- Full automated regression/integration suite passed before release preparation.
- Server validation confirmed the external database workflow and final summary schema.

## [0.1.1] - 2026

- Maintenance release for the original PlasRiP architecture.
- Bundled host-range model with external PTU sketch library.

## [0.1.0] - 2026-04-23

### Added
- Initial standardized release layout and packaging workflow.
- Release helper scripts, version synchronization, CI smoke tests, and conda packaging support.
