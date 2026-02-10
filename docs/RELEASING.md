# Releasing ZNA

This is the single reference for all ZNA publishing: PyPI, Bioconda, and GitHub Actions.

## Quick Reference

```
PyPI release:    ./scripts/release.sh X.Y.Z
Bioconda update: ./scripts/update-conda-sha.sh X.Y.Z  →  PR to bioconda-recipes
```

---

## 1. Release to PyPI

### One-time setup

Configure PyPI **trusted publishing** so GitHub Actions can publish without API tokens:

1. Go to https://pypi.org/manage/project/zna/settings/publishing/
2. Click **"Add a new publisher"**
3. Fill in:
   - **Owner:** `mkiyer`
   - **Repository:** `zna`
   - **Workflow name:** `publish.yml`
   - **Environment:** `release`
4. Click **"Add"**

That's it. No secrets or API tokens needed.

### Cutting a release

```bash
# Make sure tests pass and working directory is clean
pytest
git status

# Run the release script (updates version, tags, pushes)
./scripts/release.sh 0.2.0
```

The script will:
1. Update `__version__` in `src/zna/__init__.py`
2. Update version in `conda/meta.yaml`
3. Show a diff for your review
4. Commit, tag `v0.2.0`, and push

GitHub Actions then builds wheels for Linux/macOS/Windows × Python 3.10–3.14 and
publishes to PyPI automatically.

### Verify the release

```bash
# Check GitHub Actions
open https://github.com/mkiyer/zna/actions

# After Actions complete, verify on PyPI
open https://pypi.org/project/zna/

# Test install
pip install --upgrade zna
zna --help
```

### Manual PyPI upload (fallback)

If GitHub Actions fails:

```bash
pip install --upgrade build twine
rm -rf dist/ build/
python -m build
twine check dist/*
twine upload dist/*
# Username: __token__
# Password: your pypi- token
```

---

## 2. Update Bioconda

Do this **after** the PyPI release is live and the GitHub tag exists.

### Update the SHA256 hash

```bash
./scripts/update-conda-sha.sh 0.2.0
```

This downloads the release tarball, computes its SHA256, and patches `conda/meta.yaml`.

### Test the conda build locally (optional)

```bash
conda build conda/
```

### Commit and push

```bash
git add conda/meta.yaml
git commit -m "Update conda SHA256 for v0.2.0"
git push origin main
```

### Submit PR to bioconda-recipes

If this is the **first submission**, fork bioconda-recipes and create the recipe:

```bash
git clone https://github.com/mkiyer/bioconda-recipes.git
cd bioconda-recipes
git remote add upstream https://github.com/bioconda/bioconda-recipes.git
git fetch upstream
git checkout -b add-zna upstream/master
mkdir -p recipes/zna
cp /path/to/zna/conda/meta.yaml recipes/zna/
git add recipes/zna/
git commit -m "Add zna: high-performance nucleic acid compression"
git push origin add-zna
```

Then open a PR at https://github.com/bioconda/bioconda-recipes/pulls.

For **version updates** (including resubmitting after reviewer feedback):

```bash
cd bioconda-recipes
git fetch upstream
git checkout master && git merge upstream/master

# Create or reuse your branch
git checkout -b update-zna-0.2.0

# Copy the updated meta.yaml
cp /path/to/zna/conda/meta.yaml recipes/zna/meta.yaml

git add recipes/zna/meta.yaml
git commit -m "Update zna to 0.2.0"
git push origin update-zna-0.2.0
```

Open a new PR (or force-push to the existing branch if you're addressing reviewer
feedback on an open PR).

### Addressing reviewer feedback

When a reviewer requests changes (e.g., removing `build.sh`, inlining the build
command), make the changes in your ZNA repo, then copy the updated
`conda/meta.yaml` into your bioconda-recipes fork and push:

```bash
# In your bioconda-recipes fork, on the PR branch:
cp /path/to/zna/conda/meta.yaml recipes/zna/meta.yaml
git add recipes/zna/meta.yaml
git commit -m "Address reviewer feedback"
git push origin <your-branch>
```

The PR will update automatically.

---

## 3. Version Numbering

Follow [Semantic Versioning](https://semver.org/):

| Change type               | Example      |
|---------------------------|--------------|
| Breaking API change       | 1.0.0        |
| New feature, compatible   | 0.2.0        |
| Bug fix                   | 0.1.1        |

Version is stored in one place: `src/zna/__init__.py`. The build system
(`pyproject.toml`) and release script read it from there. `conda/meta.yaml`
has its own copy that the release script also updates.

---

## 4. Troubleshooting

### GitHub Actions workflow doesn't trigger
- The tag **must** start with `v` (e.g., `v0.2.0`, not `0.2.0`)
- Check that `.github/workflows/publish.yml` exists on the tagged commit

### PyPI publish step fails with 403 / authentication error
- Trusted publishing not configured — follow the one-time setup above
- Or the workflow `environment` name doesn't match (must be `release`)

### Need to redo a tag

```bash
git tag -d v0.2.0
git push origin :refs/tags/v0.2.0
git tag -a v0.2.0 -m "Release 0.2.0"
git push origin v0.2.0
```

### Conda build fails locally
- Ensure `cmake >= 3.15`, a C++ compiler, and `conda-build` are installed
- Check that `conda/meta.yaml` has the correct SHA256 for the release tarball

### Bioconda CI fails
- Lint your recipe: `bioconda-utils lint --git-range HEAD`
- Build with Docker: `bioconda-utils build --docker --mulled-test --git-range HEAD`

---

## Files Involved

| File | Purpose |
|------|---------|
| `src/zna/__init__.py` | Single source of truth for version |
| `pyproject.toml` | Build config; reads version dynamically |
| `conda/meta.yaml` | Conda recipe (version + SHA256) |
| `.github/workflows/publish.yml` | CI: build wheels + publish to PyPI |
| `scripts/release.sh` | Automates version bump, tag, push |
| `scripts/update-conda-sha.sh` | Updates conda SHA256 after release |

---

## Release Checklist

```
[ ] Tests pass: pytest
[ ] Working directory clean: git status
[ ] Run release: ./scripts/release.sh X.Y.Z
[ ] GitHub Actions green: check Actions tab
[ ] PyPI updated: https://pypi.org/project/zna/
[ ] Update conda SHA: ./scripts/update-conda-sha.sh X.Y.Z
[ ] Commit SHA update and push
[ ] Submit/update bioconda-recipes PR
```
