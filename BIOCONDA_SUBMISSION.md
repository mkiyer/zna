# Bioconda Submission Guide for ZNA

This guide walks through the process of submitting the ZNA package to Bioconda.

## Prerequisites

### 1. Prepare Your Repository

- [x] Source code on GitHub (public repository)
- [ ] Release tag created (e.g., `v0.1.0`)
- [ ] LICENSE file present
- [ ] README.md with clear documentation
- [ ] Working conda recipe (see `conda/meta.yaml`)

### 2. Install Required Tools

```bash
# Install conda-build and bioconda-utils
conda install -c conda-forge conda-build anaconda-client
conda install -c bioconda bioconda-utils

# Or create a dedicated environment
conda create -n bioconda-build conda-build anaconda-client bioconda-utils
conda activate bioconda-build
```

## Step-by-Step Submission Process

### Step 1: Create GitHub Release

1. **Tag your release:**
   ```bash
   git tag -a v0.1.0 -m "Release version 0.1.0"
   git push origin v0.1.0
   ```

2. **Create release on GitHub:**
   - Go to your repo → Releases → Create a new release
   - Select the v0.1.0 tag
   - Title: "ZNA v0.1.0"
   - Description: Copy from README or CHANGELOG
   - Publish release

3. **Get the SHA256 hash:**
   ```bash
   # Download the release tarball
   wget https://github.com/yourusername/zna/archive/v0.1.0.tar.gz
   
   # Calculate SHA256
   sha256sum v0.1.0.tar.gz
   # Or on macOS:
   shasum -a 256 v0.1.0.tar.gz
   ```

4. **Update meta.yaml:**
   - Replace `yourusername` with your actual GitHub username
   - Update the `sha256` field with the calculated hash

### Step 2: Test Locally

```bash
# Navigate to your project directory
cd /path/to/zna

# Test the conda build locally
conda build conda/

# If build succeeds, test the package
conda create -n test-zna python=3.14
conda activate test-zna
conda install --use-local zna

# Run tests
zna --help
python -c "from zna.core import is_accelerated; print(f'Accelerated: {is_accelerated()}')"
pytest  # If you have tests in the package

# Clean up
conda deactivate
conda env remove -n test-zna
```

### Step 3: Fork bioconda-recipes

1. **Fork the repository:**
   - Go to https://github.com/bioconda/bioconda-recipes
   - Click "Fork" button
   - Clone your fork:
     ```bash
     git clone https://github.com/yourusername/bioconda-recipes.git
     cd bioconda-recipes
     ```

2. **Add upstream remote:**
   ```bash
   git remote add upstream https://github.com/bioconda/bioconda-recipes.git
   git fetch upstream
   ```

### Step 4: Create Recipe Branch

```bash
# Create a new branch for your recipe
git checkout -b add-zna

# Create the recipe directory
mkdir -p recipes/zna

# Copy your recipe files
cp /path/to/zna/conda/meta.yaml recipes/zna/
cp /path/to/zna/conda/build.sh recipes/zna/
# Note: Do NOT include bld.bat - Bioconda doesn't build Windows packages

# Add and commit
git add recipes/zna/
git commit -m "Add zna recipe: high-performance nucleic acid compression"
```

### Step 5: Test with Bioconda Utils

```bash
# Lint the recipe (checks for common issues)
bioconda-utils lint --git-range HEAD

# Build and test (this uses Docker/containers)
bioconda-utils build --docker --mulled-test --git-range HEAD

# If you don't have Docker, test locally
conda build recipes/zna
```

### Step 6: Push and Create Pull Request

```bash
# Push your branch
git push origin add-zna
```

1. **Create PR on GitHub:**
   - Go to https://github.com/bioconda/bioconda-recipes
   - Click "Pull requests" → "New pull request"
   - Click "compare across forks"
   - Select your fork and branch
   - Title: "Add zna: high-performance nucleic acid compression"
   - Description:
     ```markdown
     ## Description
     Adds zna package - a high-performance binary format for compressed nucleic acid sequences.
     
     ## Performance
     - 135 MB/s roundtrip throughput
     - 2.8+ GB/s for long read processing
     - 3.7-4.0x compression ratio
     
     ## Features
     - C++ acceleration with pure Python fallback
     - Supports single-end, paired-end, and interleaved reads
     - Block-based streaming architecture
     
     ## Checklist
     - [x] Recipe follows Bioconda guidelines
     - [x] Builds successfully on Linux
     - [x] Includes tests
     - [x] License file included
     - [x] C++ compiler dependencies specified
     ```

2. **Submit the PR**

### Step 7: Address Review Feedback

Bioconda reviewers will check:

1. **Recipe quality:**
   - Correct dependencies
   - Proper version pinning
   - Test coverage
   - License compatibility

2. **Build tests:**
   - Linux builds pass
   - macOS builds pass (if applicable)
   - Tests execute successfully

3. **Common feedback items:**
   - Pin versions more strictly if needed
   - Add missing run dependencies
   - Improve test coverage
   - Update documentation links

**Responding to feedback:**
```bash
# Make changes locally
cd bioconda-recipes
git checkout add-zna

# Edit files as requested
vim recipes/zna/meta.yaml

# Commit and push
git add recipes/zna/
git commit -m "Address reviewer feedback: [describe changes]"
git push origin add-zna
```

### Step 8: Automatic Publication

Once approved and merged:
1. Bioconda CI builds packages for all platforms
2. Packages uploaded to anaconda.org/bioconda channel
3. Available within hours: `conda install -c bioconda zna`

## Recipe Maintenance

### Updating the Package

When you release a new version:

1. **Create new release on GitHub**
2. **Update recipe in bioconda-recipes:**
   ```bash
   cd bioconda-recipes
   git checkout master
   git pull upstream master
   git checkout -b update-zna-0.2.0
   
   # Update meta.yaml with new version and SHA256
   vim recipes/zna/meta.yaml
   
   git add recipes/zna/meta.yaml
   git commit -m "Update zna to version 0.2.0"
   git push origin update-zna-0.2.0
   ```
3. **Create PR** (same as initial submission)

### Common Issues and Solutions

**Issue: Build fails with compiler errors**
```yaml
# Ensure you have all compiler dependencies
requirements:
  build:
    - {{ compiler('c') }}
    - {{ compiler('cxx') }}
    - cmake >=3.15
```

**Issue: Import errors during tests**
```yaml
# Make sure all runtime dependencies are listed
requirements:
  run:
    - python
    - zstandard
```

**Issue: C++ extension not building**
```yaml
# Ensure host dependencies include build tools
requirements:
  host:
    - python
    - pip
    - scikit-build-core >=0.10
    - nanobind >=2.0
```

**Issue: Tests fail**
```yaml
# Verify test commands work
test:
  commands:
    - zna --help
    - python -c "import zna; from zna.core import is_accelerated"
```

## Additional Resources

- **Bioconda Guidelines**: https://bioconda.github.io/contributor/guidelines.html
- **Bioconda Docs**: https://bioconda.github.io/
- **Example Recipes**: https://github.com/bioconda/bioconda-recipes/tree/master/recipes
- **Similar C++ packages**: Look at recipes for packages with nanobind/pybind11

## Bioconda Requirements Checklist

- [ ] Package name is lowercase
- [ ] Recipe in `recipes/packagename/meta.yaml`
- [ ] Build script for Unix (`build.sh`) - NO Windows .bat files
- [ ] All dependencies specified correctly
- [ ] Tests included and pass
- [ ] License compatible with Bioconda (GPL-3.0 ✓)
- [ ] Recipe maintainer specified
- [ ] Source URL points to stable release
- [ ] SHA256 hash is correct
- [ ] Version pinning appropriate
- [ ] Works on Linux (required)
- [ ] Works on macOS (recommended)
- [ ] Documentation links valid

## Contact

If you need help:
- **Bioconda Gitter**: https://gitter.im/bioconda/Lobby
- **GitHub Issues**: https://github.com/bioconda/bioconda-recipes/issues
- **Bioconda Documentation**: https://bioconda.github.io/

Good luck with your submission! 🎉
