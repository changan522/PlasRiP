# PlasRiP tagging and release convention

## Software release

For this release:

- package version: `0.3.0`
- Git tag: `v0.3.0`
- GitHub release title: `PlasRiP v0.3.0`

Use an annotated software tag:

```bash
git tag -a v0.3.0 -m "PlasRiP v0.3.0"
```

Recommended release sequence:

```bash
python scripts/check_version_sync.py
pytest -q
git status
git add .
git commit -m "release: v0.3.0"
git tag -a v0.3.0 -m "PlasRiP v0.3.0"
git push origin main
git push origin v0.3.0
```

Recommended GitHub software-release assets:

- `plasrip-0.3.0-py3-none-any.whl`
- GitHub-generated source archives

Use `RELEASE_NOTES_v0.3.0.md` as the basis for the release description.

## Reference database

The reference database is released separately:

- database tag: `database-v3.0`
- release title: `PlasRiP Reference Database v3.0`
- release page: https://github.com/changan522/PlasRiP/releases/tag/database-v3.0

Do not bundle the MPC/PGG reference database inside the software wheel or source tree.
