#!/usr/bin/env bash
set -Eeuo pipefail
set -o errtrace
set -o functrace

echo "BASH_VERSION=${BASH_VERSION:-<unset>}"
echo "shopt -o errtrace => $(set -o | grep -E '^errtrace|^functrace')"
echo "trap -p BEFORE =>"
trap -p

# minimal ERR/DEBUG to show they fire at all
trap 'echo "[SELFTEST DEBUG] about to run: $BASH_COMMAND" >&2' DEBUG
trap 'st=$?; echo "[SELFTEST  ERR ] last: $BASH_COMMAND (exit $st)" >&2' ERR

echo "[SELFTEST] running failing command..."
false # <- simple failing command

echo "[SELFTEST] should not reach here if -e is on"
