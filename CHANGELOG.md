# Changelog

All notable changes to this project will be documented in this file.

The format is based on Keep a Changelog, and this project follows Semantic Versioning.

## [0.1.0] - 2026-04-23

### Added
- Standardized conda release layout for packaging and publication.
- Release helper scripts for build, local install, version checks, and Anaconda upload.
- A canonical `VERSION` file and a version consistency checker.
- GitHub-facing release notes and tag naming guidance for `v0.1.0`.
- A Linux CI workflow for smoke testing installation and CLI entry points.

### Changed
- Synchronized the package version across `pyproject.toml`, `src/plasrip/__init__.py`, `src/plasrip/cli.py`, and `conda.recipe/meta.yaml`.
- Expanded release documentation for GitHub and Anaconda.org publication.
- Cleaned the release tree by removing caches and non-release build outputs.

### Fixed
- Preserved a stable fallback version string in the CLI for source-based execution.
- Kept conda build and smoke-test commands aligned with the release directory structure.
