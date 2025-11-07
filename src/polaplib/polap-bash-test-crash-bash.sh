#!/usr/bin/env bash
# Version: v0.1.0
# Demo: Bash crash with FAILED CMD + call stack.
set -Eeuo pipefail

# expected env: _POLAPLIB_DIR points to your libs; fallback to script dir
_POLAPLIB_DIR="${_POLAPLIB_DIR:-"$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"}"
source "${_POLAPLIB_DIR}/polap-lib-logcallstack.sh"
polap_trap_enable

work="${1:-demo-bash}"
mkdir -p "$work/A" "$work/B"
: >"$work/A/file.txt"

# First link (ok or preexisting)
ln -sf ../A/file.txt "$work/B/link.txt" || true
# Second link without -f will fail with "File exists" → triggers ERR
ln -s ../A/file.txt "$work/B/link.txt"
