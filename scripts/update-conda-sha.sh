#!/bin/bash
# Update (or just verify) the SHA256 in conda/meta.yaml after a GitHub release.
#
# Usage: ./scripts/update-conda-sha.sh <version> [--check]
#   ./scripts/update-conda-sha.sh 0.4.0            # download, hash, rewrite meta.yaml
#   ./scripts/update-conda-sha.sh 0.4.0 --check    # verify only; exits 1 on mismatch
#
# --check exists because nothing else in the release verifies that conda/meta.yaml's
# sha256 actually matches the tarball at its own URL. Until it was added, bioconda CI
# was the first thing to notice -- which is how v0.4.0 shipped a recipe whose hash had
# gone stale after the tag was re-cut.

set -euo pipefail

if [ -z "${1:-}" ]; then
    echo "Error: Version number required"
    echo "Usage: ./scripts/update-conda-sha.sh <version> [--check]"
    exit 1
fi

VERSION=$1
MODE=${2:-update}
META="conda/meta.yaml"

if [ ! -f "$META" ]; then
    echo "❌ $META not found. Run this from the repository root."
    exit 1
fi

USERNAME=$(grep 'set username' "$META" | sed 's/.*"\(.*\)".*/\1/')
PACKAGE="zna"
URL="https://github.com/${USERNAME}/${PACKAGE}/archive/v${VERSION}.tar.gz"

echo "📥 Downloading release tarball..."
echo "URL: ${URL}"

TMP_FILE=$(mktemp)
trap 'rm -f "$TMP_FILE"' EXIT

# -f is load-bearing. Without it a missing tag returns HTTP 404 with a 14-byte
# "404: Not Found" BODY, curl exits 0, and this script hashes the error page and
# reports success -- writing d5558cd419c8d4... into the recipe. Verified.
if ! curl -fsSL --retry 3 "${URL}" -o "${TMP_FILE}"; then
    echo ""
    echo "❌ Download failed. Does the tag v${VERSION} exist and is it pushed?"
    echo "   Check: git ls-remote --tags origin v${VERSION}"
    exit 1
fi

# Belt and braces: a proxy or a rate-limit page can still return HTTP 200 with HTML.
# A release tarball starts with the gzip magic bytes 1f 8b.
if [ "$(head -c 2 "${TMP_FILE}" | xxd -p)" != "1f8b" ]; then
    echo ""
    echo "❌ Downloaded $(wc -c < "${TMP_FILE}" | tr -d ' ') bytes that are not a gzip archive:"
    head -c 200 "${TMP_FILE}"
    echo ""
    exit 1
fi

SHA256=$(shasum -a 256 "${TMP_FILE}" | cut -d' ' -f1)
CURRENT=$(grep -E '^[[:space:]]*sha256:' "$META" | head -1 | awk '{print $2}')

echo ""
echo "✅ SHA256: ${SHA256}"
echo ""

if [ "$MODE" = "--check" ]; then
    if [ "$SHA256" = "$CURRENT" ]; then
        echo "✅ ${META} matches the tarball at v${VERSION}."
        exit 0
    fi
    echo "❌ ${META} is STALE."
    echo "   recipe says: ${CURRENT}"
    echo "   tarball is:  ${SHA256}"
    echo ""
    echo "   Re-run without --check to fix. If the tag was deleted and re-cut, the"
    echo "   tarball bytes changed and every copy of this recipe must be re-hashed --"
    echo "   including the one in bioconda-recipes."
    exit 1
fi

# Exactly one sha256 line must exist, because the rewrite below is a blanket
# substitution. If the recipe ever grows a second source, this is the guard that
# stops both being overwritten with the same hash.
N=$(grep -cE '^[[:space:]]*sha256:' "$META")
if [ "$N" -ne 1 ]; then
    echo "❌ Expected exactly one sha256: line in ${META}, found ${N}."
    echo "   Refusing to rewrite; update it by hand."
    exit 1
fi

if [ "$SHA256" = "$CURRENT" ]; then
    echo "ℹ️  ${META} is already correct; nothing to change."
    exit 0
fi

echo "📝 Updating ${META}..."
echo "   ${CURRENT}"
echo "-> ${SHA256}"
# POSIX class, not \s: BSD sed (the macOS default) has no \s in ERE and treats it as a
# literal 's', so this pattern silently matched nothing on the maintainer's own machine
# while working fine under GNU sed on Linux. BSD *grep* does accept \s, which is why the
# reads below looked correct and only the rewrite failed. The verify block caught it.
sed -i.bak -E "s|^([[:space:]]*)sha256:.*|\1sha256: ${SHA256}|" "$META"
rm -f "${META}.bak"

# Confirm the edit actually landed, rather than trusting sed.
WROTE=$(grep -E '^[[:space:]]*sha256:' "$META" | head -1 | awk '{print $2}')
if [ "$WROTE" != "$SHA256" ]; then
    echo "❌ Rewrite did not take: ${META} still says ${WROTE}"
    exit 1
fi

echo ""
echo "✅ Updated ${META}"
echo ""
echo "📋 Next steps:"
echo "   1. Review:  git diff ${META}"
echo "   2. Commit:  git add ${META} && git commit -m 'Update conda SHA256 for v${VERSION}'"
echo "   3. Push:    git push origin main"
echo "   4. Bioconda: the autobump bot opens a version-only PR that keeps the PREVIOUS"
echo "      release's recipe body. Overwrite recipes/zna/meta.yaml with THIS file --"
echo "      do not merely re-hash the bot's branch, or the new recipe's tests and"
echo "      dependencies never ship. See docs/RELEASING.md section 2."
