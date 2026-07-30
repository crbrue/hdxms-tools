#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="${1:-$(pwd)}"
SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/resources/styles"
DEST_DIR="${REPO_ROOT}/resources/styles"

if [[ ! -d "${REPO_ROOT}/src/hdxms" ]]; then
    echo "ERROR: ${REPO_ROOT} does not look like the hdx-ms-tools repository."
    echo "Usage: bash install_madp_hdx_palette.sh /path/to/hdx-ms-tools"
    exit 1
fi

mkdir -p "${DEST_DIR}"
cp -v "${SOURCE_DIR}/"* "${DEST_DIR}/"

echo
echo "Installed MADP HDX Palette v1.0 into:"
echo "  ${DEST_DIR}"
echo
echo "Canonical source:"
echo "  ${DEST_DIR}/hdx_palette.yaml"