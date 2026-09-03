PlasRiP conda recipe notes (v0.2.0)

- Host-range model assets under `src/plasrip/data/hostrange/` remain packaged with the software.
- CONJScan assets remain packaged for mobility analysis.
- PTU/PGG reference metadata and Mash sketches are NOT bundled with the conda package.
- PTU-related workflows require the external PlasRiP Reference Database v2.0 via `--database` or `PLASRIP_DATABASE`.
- The external database is versioned independently from the software package.
