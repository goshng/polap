#!/usr/bin/env bash
# File: polap-bash-run-command-backend-smoketest.sh
# Version: v0.2.4
# Purpose: Smoke test for POLAP_RUN_BACKEND selection.

set -Eeuo pipefail
LIBDIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
source "${LIBDIR}/polap-lib-run-command.sh"

TMP=$(mktemp -d)
export OUTDIR="${TMP}/out"
polap__ensure_dir "$OUTDIR"

echo "[smoke] OUTDIR=$OUTDIR"

# auto prefers argv, then simple, then direct
POLAP_RUN_BACKEND=auto polap_run_wrapper "auto_ok" /bin/true

# explicit argv
POLAP_RUN_BACKEND=argv polap_run_wrapper "argv_ok" echo "argv backend ok"

# explicit simple: requires a shell string with a pipe
POLAP_RUN_BACKEND=simple polap_run_wrapper "simple_ok" "printf 'ok\n' | wc -l"

# explicit direct: runs argv without shell
POLAP_RUN_BACKEND=direct polap_run_wrapper "direct_ok" bash -c 'echo direct'

echo "[smoke] done. Check:"
find "$OUTDIR/.polap-run" -maxdepth 2 -name .ok -printf '%P\n' | sort
