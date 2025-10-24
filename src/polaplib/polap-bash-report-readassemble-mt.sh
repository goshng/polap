#!/usr/bin/env bash
# polap-bash-report-readassemble-mt.sh
# Version: v0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Build a JSON-only report for MT pipeline (panel -> PT-gating -> allvsall -> busco -> miniasm -> mt{0,1,2,3}).
# Input:  --base-dir points to .../polap-readassemble/mtseed
# Output: <base-dir>/report/mt-report.json
#
set -euo pipefail
IFS=$'\n\t'

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --base-dir /path/to/polap-readassemble/mtseed
EOF
}

BASE_DIR=""

while (($#)); do
	case "$1" in
	--base-dir)
		BASE_DIR="${2:?}"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERR] Unknown option: $1" >&2
		usage
		exit 2
		;;
	esac
done

[[ -n "$BASE_DIR" && -d "$BASE_DIR" ]] || {
	echo "[ERR] --base-dir invalid: $BASE_DIR" >&2
	exit 2
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POLAPLIB_DIR="${POLAPLIB_DIR:-$SCRIPT_DIR}"
PY_TOOL="${POLAPLIB_DIR}/scripts/report_readassemble_mt.py"

[[ -f "$PY_TOOL" ]] || {
	echo "[ERR] Missing Python helper: $PY_TOOL" >&2
	exit 3
}

OUT_DIR="${BASE_DIR}/report"
mkdir -p "$OUT_DIR"
OUT_JSON="${OUT_DIR}/mt-report.json"

python3 "$PY_TOOL" --base-dir "$BASE_DIR" --out "$OUT_JSON"

echo "[OK] MT JSON: $OUT_JSON"
