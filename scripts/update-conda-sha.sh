#!/bin/bash
# Update SHA256 in conda/meta.yaml after GitHub release
# Usage: ./scripts/update-conda-sha.sh <version>
# Example: ./scripts/update-conda-sha.sh 0.2.0

set -e

if [ -z "$1" ]; then
    echo "Error: Version number required"
    echo "Usage: ./scripts/update-conda-sha.sh <version>"
    exit 1
fi

VERSION=$1
USERNAME=$(grep 'set username' conda/meta.yaml | sed 's/.*"\(.*\)".*/\1/')
PACKAGE="zna"
URL="https://github.com/${USERNAME}/${PACKAGE}/archive/v${VERSION}.tar.gz"

echo "📥 Downloading release tarball..."
echo "URL: ${URL}"

# Download and calculate SHA256
TMP_FILE=$(mktemp)
curl -L "${URL}" -o "${TMP_FILE}"
SHA256=$(shasum -a 256 "${TMP_FILE}" | cut -d' ' -f1)
rm "${TMP_FILE}"

echo ""
echo "✅ SHA256: ${SHA256}"
echo ""

# Update meta.yaml
echo "📝 Updating conda/meta.yaml..."
sed -i.bak "s/sha256: .*/sha256: ${SHA256}/" conda/meta.yaml
rm conda/meta.yaml.bak

echo ""
echo "✅ Updated conda/meta.yaml"
echo ""
echo "📋 Next steps:"
echo "   1. Review changes: git diff conda/meta.yaml"
echo "   2. Test build: conda build conda/"
echo "   3. Commit: git add conda/meta.yaml && git commit -m 'Update conda SHA256 for v${VERSION}'"
echo "   4. Push: git push origin main"
echo "   5. Submit PR to bioconda-recipes repository"
