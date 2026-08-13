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

**The release must be cut from `main`.** `release.sh` pushes `main` and then tags
`HEAD`, so running it from a feature branch would push whatever local `main` happens to
point at and then tag a commit `main` does not contain — a release whose tag and whose
branch disagree. The script now refuses rather than doing that, but the merge is still
yours to do:

```bash
# 1. Verify on the release branch, in BOTH environments (see below)
pytest

# 2. Merge to main
git checkout main
git merge zna-0.4.0-hardening
git status                      # must be clean

# 3. Cut it
./scripts/release.sh 0.4.0
```

The script will:
1. Refuse unless the tree is clean and you are on `main`
2. Update `__version__` in `src/zna/__init__.py` and the version in `conda/meta.yaml` —
   or say so and carry on if they are already at the target, which is the normal case
   when development has bumped `__init__.py` ahead of the tag
3. Show a diff for your review, and prompt before doing anything irreversible
4. Commit (if there is anything to commit), tag `v0.4.0`, and push both

#### Check CI is green on the branch before you tag

`.github/workflows/ci.yml` builds and tests on Linux, macOS and Windows on every push.
**Wait for it.** A local test run proves nothing about MSVC: 0.4.0's first tag had to be
deleted because the Windows wheel failed on a compiler builtin that a macOS developer
machine never even compiles. The tag is the expensive place to discover that.

#### Verify in both environments first

The compiled merge extension is the thing most likely to be silently absent, and it
fails *quietly correct* — the reference kernel is ~50x slower and gives identical
answers, so a broken build looks like a slow node rather than a mistake.

```bash
# compiled: the configuration that ships
/path/to/envs/zna_merge/bin/python -m pytest        # expect 0 skips in the merge suite

# extension-less: proves the reference oracle is not silently untested
/path/to/envs/zna/bin/python -m pytest              # expect ~56 skips
```

A rebuild is `pip install -e . --no-build-isolation`, and **grep its output for
`Successfully built`** — a CMake failure ends with an unremarkable pip banner, leaves
the old `.so` installed, and lets the whole suite run green against stale code.

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
[ ] CI green on the branch, all three platforms  <-- MSVC is only checked here
[ ] Tests pass in the COMPILED env (0 merge skips)
[ ] Tests pass in the EXTENSION-LESS env (~56 skips, oracle exercised)
[ ] Rebuild output contained "Successfully built"
[ ] Merged to main:  git checkout main && git merge <branch>
[ ] Working directory clean: git status
[ ] Run release: ./scripts/release.sh X.Y.Z
[ ] GitHub Actions green: check Actions tab
[ ] PyPI updated: https://pypi.org/project/zna/
[ ] Wheel smoke test: python scripts/wheel_smoketest.py  (see below)
[ ] Update conda SHA: ./scripts/update-conda-sha.sh X.Y.Z
[ ] Commit SHA update and push
[ ] Submit/update bioconda-recipes PR
```

### Why the wheel needs its own check

ZNA builds **two** extensions — `zna._accel` (the codec) and `zna.merge._accel` (the
overlap scan). They are different CMake targets that both produce a file named
`_accel.cpython-*.so`, and they have been installed into each other's place before; when
that happened `zna.is_accelerated()` returned False while every test still passed. A
wheel is the one artifact where nobody would notice, so check the installed package
rather than the source tree:

```bash
pip install --upgrade zna
python -c "
import zna
from zna.merge.backend import available_merge_backends
assert zna.is_accelerated(), 'codec extension missing from the wheel'
assert 'accel' in available_merge_backends(), 'merge extension missing from the wheel'
print('ok', zna.__version__, available_merge_backends())
"
```
