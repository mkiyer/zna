# Quick PyPI Publishing Guide for ZNA

## One-Time Setup

1. **Create PyPI account**: https://pypi.org/account/register/
2. **Generate API token**: https://pypi.org/manage/account/ → "API tokens" → "Add API token"
3. **Save token** (starts with `pypi-...`)

## Publishing Steps

### 1. Clean and Build

```bash
cd /Users/mkiyer/proj/zna
rm -rf dist/ build/
python -m build
```

### 2. Check Distribution

```bash
twine check dist/*
```

### 3. Upload to PyPI

```bash
twine upload dist/*
```

When prompted:
- Username: `__token__`
- Password: `pypi-YourActualTokenHere` (paste your token)

### 4. Verify

```bash
# Wait ~1 minute for PyPI to process
pip install --upgrade zna
zna --help
```

## For Future Releases

1. Update version in `pyproject.toml`
2. Commit changes:
   ```bash
   git add pyproject.toml
   git commit -m "Bump version to X.Y.Z"
   git tag vX.Y.Z
   git push origin main --tags
   ```
3. Rebuild and upload:
   ```bash
   rm -rf dist/ build/
   python -m build
   twine check dist/*
   twine upload dist/*
   ```

## Your Current Build

✅ `dist/zna-0.1.2.tar.gz` - Source distribution
✅ `dist/zna-0.1.2-cp314-cp314-macosx_15_0_arm64.whl` - macOS ARM64 wheel

Both passed validation and are ready for upload!

**Note on file extensions:** ZNA uses `.zna` for all files (both compressed and uncompressed). Compression is enabled by default (Zstd level 3) and can be verified with `zna inspect`.

---

**Ready to publish?** Run: `twine upload dist/*`
