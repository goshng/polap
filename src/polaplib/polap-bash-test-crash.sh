#!/usr/bin/env bash
# polap-bash-test-crash.sh
# Version: v0.1.0
set -Eeuo pipefail
: "${_POLAPLIB_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_POLAPLIB_DIR}/polap-lib-logcallstack.sh"
polap_trap_enable

work="${_POLAPLIB_DIR}/_demo_bash"
mkdir -p "$work/A" "$work/B"
: >"$work/A/file.txt"

# First link (fine or overwritten)
ln -sf ../A/file.txt "$work/B/link.txt" || true
# Second link without -f will fail with "File exists" ⇒ triggers ERR/stack
ln -s ../A/file.txt "$work/B/link.txt"
