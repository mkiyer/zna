# GitHub Actions Setup for Automated Releases

To enable automated PyPI publishing, you need to configure a GitHub secret.

## Step 1: Create PyPI API Token

1. Go to https://pypi.org/manage/account/token/
2. Click "Add API token"
3. Token name: `zna-github-actions`
4. Scope: "Entire account" (or just "Project: zna" if available)
5. Click "Add token"
6. **Copy the token** (starts with `pypi-...`)
   - ⚠️ You can only see this once!

## Step 2: Add Secret to GitHub

1. Go to your repository: https://github.com/mkiyer/zna
2. Click **Settings** → **Secrets and variables** → **Actions**
3. Click **"New repository secret"**
4. Set:
   - Name: `PYPI_TOKEN`
   - Value: Paste your PyPI token (the `pypi-...` string)
5. Click **"Add secret"**

## Step 3: Verify Setup

To test the automation:

```bash
# Create a test release
./scripts/release.sh 0.1.4

# Then check GitHub Actions:
# https://github.com/mkiyer/zna/actions
```

The workflow will:
1. Trigger on the `v0.1.4` tag
2. Build wheels for all platforms
3. Publish to PyPI automatically

## Troubleshooting

### "Secret PYPI_TOKEN not found"
- Verify the secret name is exactly `PYPI_TOKEN` (case-sensitive)
- Ensure you're adding it as a **repository secret**, not an environment secret

### "Invalid or expired token"
- Generate a new token on PyPI
- Update the GitHub secret with the new token

### Workflow doesn't trigger
- Ensure the tag starts with `v` (e.g., `v0.1.4`, not `0.1.4`)
- Check the workflow file: `.github/workflows/publish.yml`

## Security Notes

- ✅ API tokens are more secure than username/password
- ✅ Tokens can be scoped to specific projects
- ✅ Tokens can be revoked at any time
- ✅ GitHub Secrets are encrypted and never exposed in logs

## Manual Publishing (Fallback)

If GitHub Actions is unavailable, you can still publish manually:

```bash
python -m build
twine upload dist/*
```
