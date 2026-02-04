# Publishing ZNA to PyPI

> **⚡ Quick Start:** For streamlined releases, see [RELEASING.md](./RELEASING.md)
>
> This document describes the **manual publishing process**. For most releases, use the automated workflow instead.

## Automated Release (Recommended)

Use the automated release script:

```bash
./scripts/release.sh 0.2.0
```

This automatically:
- Updates version in all files
- Creates git tag
- Triggers GitHub Actions to build and publish to PyPI

See [RELEASING.md](./RELEASING.md) for details.

---

## Manual Publishing (Fallback)

This guide explains how to manually build and publish if the automated process is unavailable.

## Prerequisites

1. **Create PyPI Account**
   - Register at https://pypi.org/account/register/
   - Verify your email address
   - Enable 2FA (required for new projects)

2. **Generate API Token**
   - Go to https://pypi.org/manage/account/
   - Scroll to "API tokens"
   - Click "Add API token"
   - Set scope to "Entire account" (or specific to zna after first upload)
   - Copy the token (starts with `pypi-`)

3. **Install Build Tools**
   ```bash
   pip install --upgrade build twine
   ```

## Building the Package

1. **Clean Previous Builds**
   ```bash
   rm -rf build/ dist/ *.egg-info
   ```

2. **Build Source Distribution and Wheel**
   ```bash
   python -m build
   ```

   This creates:
   - `dist/zna-0.1.2.tar.gz` (source distribution)
   - `dist/zna-0.1.2-*.whl` (wheel with C++ extension)

3. **Check the Distribution**
   ```bash
   twine check dist/*
   ```

## Testing on TestPyPI (Recommended First Time)

1. **Upload to TestPyPI**
   ```bash
   twine upload --repository testpypi dist/*
   ```

   When prompted, use:
   - Username: `__token__`
   - Password: Your TestPyPI token

2. **Test Installation**
   ```bash
   pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ zna
   ```

3. **Verify**
   ```bash
   zna --help
   python -c "from zna.core import is_accelerated; print(f'Accelerated: {is_accelerated()}')"
   ```

## Publishing to PyPI (Production)

1. **Upload to PyPI**
   ```bash
   twine upload dist/*
   ```

   When prompted, use:
   - Username: `__token__`
   - Password: Your PyPI token (starts with `pypi-`)

2. **Verify on PyPI**
   - Visit https://pypi.org/project/zna/
   - Check metadata, description, and links

3. **Test Installation**
   ```bash
   pip install zna
   zna --help
   ```

## Using API Token via Configuration

Create `~/.pypirc`:

```ini
[distutils]
index-servers =
    pypi
    testpypi

[pypi]
username = __token__
password = pypi-YourActualTokenHere

[testpypi]
repository = https://test.pypi.org/legacy/
username = __token__
password = pypi-YourTestPyPITokenHere
```

**Security Note:** Keep your `.pypirc` file private (permissions 600).

```bash
chmod 600 ~/.pypirc
```

Then upload without prompts:
```bash
twine upload dist/*
```

## Automated Release Workflow

For future releases:

1. **Update Version**
   - Edit version in `pyproject.toml`
   - Update `CHANGELOG.md` (if exists)

2. **Commit and Tag**
   ```bash
   git add pyproject.toml
   git commit -m "Bump version to 0.1.3"
   git tag v0.1.3
   git push origin main --tags
   ```

3. **Build and Publish**
   ```bash
   rm -rf dist/ build/
   python -m build
   twine check dist/*
   twine upload dist/*
   ```

## GitHub Actions (Optional)

Consider setting up GitHub Actions for automated publishing on release tags.

Create `.github/workflows/publish.yml`:

```yaml
name: Publish to PyPI

on:
  release:
    types: [published]

jobs:
  deploy:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v3
    - name: Set up Python
      uses: actions/setup-python@v4
      with:
        python-version: '3.10'
    - name: Install dependencies
      run: |
        pip install build twine
    - name: Build package
      run: python -m build
    - name: Publish to PyPI
      env:
        TWINE_USERNAME: __token__
        TWINE_PASSWORD: ${{ secrets.PYPI_API_TOKEN }}
      run: twine upload dist/*
```

Add your PyPI token to GitHub repository secrets as `PYPI_API_TOKEN`.

## Troubleshooting

### Build Failures

If the C++ extension fails to build:
- Ensure CMake ≥3.15 is installed
- Ensure a C++ compiler is available
- Check that nanobind is properly installed

### Upload Errors

**403 Forbidden**: Invalid token or insufficient permissions
- Verify token is correct
- Check token scope includes the package

**400 Bad Request**: Filename already exists
- Version already uploaded (increment version number)

**Package name already taken**: Choose a different name or request access

## Version Numbering

Follow [Semantic Versioning](https://semver.org/):
- **MAJOR**: Incompatible API changes
- **MINOR**: New functionality, backwards-compatible
- **PATCH**: Bug fixes, backwards-compatible

Current: `0.1.2` (Beta release)
