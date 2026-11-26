#!/usr/bin/env bash
# Writes HASHES.sha256: SHA-256 of every file of the repository.
# Verify with:  shasum -a 256 -c HASHES.sha256 (macOS)  or  sha256sum -c HASHES.sha256 (Linux)
set -euo pipefail
cd "$(dirname "$0")"
if command -v sha256sum >/dev/null; then H="sha256sum"; else H="shasum -a 256"; fi
find . -type f ! -path "./.git/*" ! -name "HASHES.sha256" ! -name ".DS_Store" ! -path "*/__pycache__/*" \
  ! -path "./.pytest_cache/*" ! -path "./figures/*" -print0 | LC_ALL=C sort -z | xargs -0 $H > HASHES.sha256
echo "HASHES.sha256: $(wc -l < HASHES.sha256) files"
