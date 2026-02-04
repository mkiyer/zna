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

# Commit changes
echo "💾 Committing version bump..."
git add src/zna/__init__.py conda/meta.yaml
git commit -m "Bump version to ${VERSION}"

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
