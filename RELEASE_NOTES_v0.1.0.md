# PlasRip v0.1.0

This is the first packaged release of **PlasRip** prepared for GitHub and Anaconda.org publication.

## Highlights
- Standard conda recipe for Linux/WSL environments
- Packaged CLI entry point: `plasrip`
- Integrated host-range prediction and PTU assignment workflow
- Bundled release scripts for build, validation, and upload

## Release contents
- synchronized version metadata (`0.1.0`)
- conda recipe under `conda.recipe/`
- release helper scripts under `scripts/`
- tests for import, CLI, and installation smoke checks

## Notes
- Windows package builds remain disabled in the conda recipe.
- External tools such as `mash`, `abricate`, `macsyfinder`, and `prodigal` are expected from conda dependencies.
- Host-range model files are not embedded in this source release and should be supplied with `--model-dir`.
