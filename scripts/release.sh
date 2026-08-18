#!/bin/bash
# Release automation script for ZNA
# Usage: ./scripts/release.sh <version>
# Example: ./scripts/release.sh 0.2.0

set -e

if [ -z "$1" ]; then
    echo "Error: Version number required"
    echo "Usage: ./scripts/release.sh <version>"
    echo "Example: ./scripts/release.sh 0.2.0"
    exit 1
fi

VERSION=$1
TAG="v${VERSION}"

echo "🚀 Releasing ZNA version ${VERSION}"
echo ""

# Check if working directory is clean
if [ -n "$(git status --porcelain)" ]; then
    echo "❌ Error: Working directory is not clean"
    echo "Please commit or stash your changes first"
    exit 1
fi

# Must be on main. This script pushes `main` and then the tag, so running it from a
# feature branch pushes whatever local main happens to point at and then tags a commit
# main does not contain -- a release whose tag and whose branch disagree.
BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [ "$BRANCH" != "main" ]; then
    echo "❌ Error: on branch '${BRANCH}', not main"
    echo "This script pushes main and tags HEAD. Merge your branch first:"
    echo "    git checkout main && git merge ${BRANCH} && ./scripts/release.sh ${VERSION}"
    exit 1
fi

# The CHANGELOG must have a stamped section for this release BEFORE tagging.
# 0.5.0 shipped with all its content still under "[Unreleased]" because nothing
# checked; this is the check.
if ! grep -q "^## \[${VERSION}\]" CHANGELOG.md; then
    echo "❌ Error: CHANGELOG.md has no '## [${VERSION}]' section."
    echo "Stamp the [Unreleased] content as [${VERSION}] - $(date +%Y-%m-%d) first."
    exit 1
fi

# Update version in __init__.py
echo "📝 Updating version in src/zna/__init__.py..."
sed -i.bak "s/__version__ = \".*\"/__version__ = \"${VERSION}\"/" src/zna/__init__.py
rm src/zna/__init__.py.bak

# Update version in conda/meta.yaml
echo "📝 Updating version in conda/meta.yaml..."
sed -i.bak "s/{% set version = \".*\" %}/{% set version = \"${VERSION}\" %}/" conda/meta.yaml
rm conda/meta.yaml.bak

# Show changes
echo ""
echo "📋 Changes:"
git diff src/zna/__init__.py conda/meta.yaml

# Confirm
echo ""
read -p "Continue with release? (y/N) " -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "❌ Release cancelled"
    git checkout src/zna/__init__.py conda/meta.yaml
    exit 1
fi

# Commit changes.
#
# The version files may ALREADY be at ${VERSION} -- development commonly bumps
# __init__.py long before the tag. `git commit` with nothing staged exits non-zero, and
# under `set -e` that aborted the release here, after the confirmation prompt and before
# the tag. Nothing to commit is a normal state, not a failure.
git add src/zna/__init__.py conda/meta.yaml
if git diff --cached --quiet; then
    echo "💾 Version already ${VERSION}; nothing to commit."
else
    echo "💾 Committing version bump..."
    git commit -m "Bump version to ${VERSION}"
fi

# Create and push tag
echo "🏷️  Creating tag ${TAG}..."
git tag -a "${TAG}" -m "Release ${VERSION}"

echo "⬆️  Pushing changes and tag..."
git push origin main
git push origin "${TAG}"

echo ""
echo "✅ Release ${VERSION} initiated!"
echo ""
echo "📦 GitHub Actions will automatically:"
echo "   1. Build wheels for all platforms"
echo "   2. Publish to PyPI"
echo ""
echo "🐍 Manual steps remaining:"
echo "   1. Update conda SHA256: ./scripts/update-conda-sha.sh ${VERSION}"
echo "   2. Test conda build locally: conda build conda/"
echo "   3. Submit PR to bioconda-recipes"
echo ""
echo "🔗 Track progress: https://github.com/$(git remote get-url origin | sed 's/.*github.com[:/]\(.*\).git/\1/')/actions"
