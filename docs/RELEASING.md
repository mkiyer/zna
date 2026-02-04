# 🚀 Quick Release Guide

Release a new version in **one command**:

```bash
./scripts/release.sh 0.2.0
```

That's it! 🎉

## What happens automatically:

1. ✅ Version updated in `__init__.py` and `meta.yaml`
2. ✅ Changes committed and tagged
3. ✅ Pushed to GitHub
4. ✅ **PyPI publishing triggered automatically**
5. ✅ Wheels built for all platforms

## After the release:

Update conda SHA256 and submit to Bioconda:

```bash
# Update SHA256
./scripts/update-conda-sha.sh 0.2.0

# Test locally
conda build conda/

# Commit and push
git add conda/meta.yaml
git commit -m "Update conda SHA256 for v0.2.0"
git push

# Submit PR to bioconda-recipes repository
```

## First-time setup:

1. **Set up PyPI token in GitHub** (one-time only)
   - See [GITHUB_ACTIONS_SETUP.md](./GITHUB_ACTIONS_SETUP.md)

2. **Make scripts executable** (already done)
   ```bash
   chmod +x scripts/*.sh
   ```

## See also:

- [RELEASE.md](./RELEASE.md) - Full documentation
- [GITHUB_ACTIONS_SETUP.md](./GITHUB_ACTIONS_SETUP.md) - GitHub Actions setup

---

**Before you release:** Ensure tests pass and working directory is clean!
