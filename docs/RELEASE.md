# Release Process

This document describes the streamlined release process for ZNA.

## Overview

The release process is now **highly automated**:
- ✅ Version stored in **one place** (`src/zna/__init__.py`)
- ✅ **One command** to release: `./scripts/release.sh <version>`
- ✅ **Automatic PyPI publishing** via GitHub Actions
- ✅ **Script-assisted conda updates**

## Prerequisites

1. **PyPI API Token** stored in GitHub Secrets as `PYPI_TOKEN`
   - Go to: Repository Settings → Secrets → Actions → New repository secret
   - Name: `PYPI_TOKEN`
   - Value: Your PyPI API token

2. **Clean working directory** (no uncommitted changes)

3. **Main branch** up to date

## Quick Release (90% Automated)

```bash
# Release version 0.2.0
./scripts/release.sh 0.2.0
```

This single command will:
1. ✅ Update version in `src/zna/__init__.py`
2. ✅ Update version in `conda/meta.yaml`
3. ✅ Show you the changes for review
4. ✅ Commit the version bump
5. ✅ Create git tag (e.g., `v0.2.0`)
6. ✅ Push to GitHub
7. ✅ **Automatically trigger PyPI publishing** (via GitHub Actions)

GitHub Actions will then:
- Build wheels for Linux, macOS, Windows
- Build for Python 3.10, 3.11, 3.12
- Publish everything to PyPI

## Conda/Bioconda Update (Semi-Automated)

After the GitHub release is created:

```bash
# Update SHA256 hash for conda
./scripts/update-conda-sha.sh 0.2.0

# Test conda build locally
conda build conda/

# Commit the SHA256 update
git add conda/meta.yaml
git commit -m "Update conda SHA256 for v0.2.0"
git push origin main
```

Then submit PR to bioconda-recipes:
1. Fork: https://github.com/bioconda/bioconda-recipes
2. Update `recipes/zna/meta.yaml` with new version and SHA256
3. Submit PR

## Manual Release (If Needed)

If you need to release manually without the script:

```bash
# 1. Update version in __init__.py
vim src/zna/__init__.py  # Change __version__ = "0.2.0"

# 2. Update conda meta.yaml
vim conda/meta.yaml  # Change {% set version = "0.2.0" %}

# 3. Commit and tag
git add src/zna/__init__.py conda/meta.yaml
git commit -m "Bump version to 0.2.0"
git tag -a v0.2.0 -m "Release 0.2.0"
git push origin main v0.2.0

# 4. GitHub Actions will auto-publish to PyPI
# 5. Follow conda steps above
```

## Version Numbering

Follow [Semantic Versioning](https://semver.org/):
- **Major** (1.0.0): Breaking changes
- **Minor** (0.2.0): New features, backward compatible
- **Patch** (0.1.1): Bug fixes

## Troubleshooting

### GitHub Actions Failed

Check workflow logs:
```bash
# Open in browser
gh workflow view publish

# Or check online
# https://github.com/mkiyer/zna/actions
```

### Need to Retag

```bash
# Delete local and remote tag
git tag -d v0.2.0
git push origin :refs/tags/v0.2.0

# Create new tag
git tag -a v0.2.0 -m "Release 0.2.0"
git push origin v0.2.0
```

### Manual PyPI Upload

If GitHub Actions fails:
```bash
python -m build
twine upload dist/*
```

## Release Checklist

- [ ] Tests passing: `pytest`
- [ ] Clean working directory: `git status`
- [ ] Run release script: `./scripts/release.sh X.Y.Z`
- [ ] Verify GitHub Actions completed: Check Actions tab
- [ ] Verify PyPI upload: https://pypi.org/project/zna/
- [ ] Update conda SHA: `./scripts/update-conda-sha.sh X.Y.Z`
- [ ] Test conda build: `conda build conda/`
- [ ] Submit bioconda PR
- [ ] Update CHANGELOG.md with release notes

## Files Involved

- `src/zna/__init__.py` - Single source of truth for version
- `pyproject.toml` - Dynamic version from __init__.py
- `conda/meta.yaml` - Version for conda builds
- `.github/workflows/publish.yml` - Automated PyPI publishing
- `scripts/release.sh` - Main release automation
- `scripts/update-conda-sha.sh` - Helper for conda updates
