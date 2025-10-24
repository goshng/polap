#!/usr/bin/env bash
# polap-bash-report-readassemble-pt.sh
# Version: v0.2.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Purpose:
#   Report KPIs for PT-1 (plastid read selection & assembly) as JSON only.
#   You provide --base-dir pointing to polap-readassemble/annotate-read-pt/.
#
# Output:
#   <base-dir>/report/pt1-report.json
#
# Usage:
#   polap-bash-report-readassemble-pt.sh --base-dir <path/to/polap-readassemble/annotate-read-pt>
#
set -euo pipefail
IFS=$'\n\t'

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --base-dir DIR

DIR should be the full path to: polap-readassemble/annotate-read-pt
Example:
  $(basename "$0") --base-dir Anthoceros_agrestis/v6/0/polap-readassemble/annotate-read-pt
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

# Locate this script's directory and the Python helper
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POLAPLIB_DIR="${POLAPLIB_DIR:-$SCRIPT_DIR}"
PY_TOOL="${POLAPLIB_DIR}/scripts/report_readassemble_pt.py"

[[ -x "$PY_TOOL" || -f "$PY_TOOL" ]] || {
	echo "[ERR] Missing Python tool: $PY_TOOL" >&2
	exit 3
}

OUT_DIR="${BASE_DIR}/report"
mkdir -p "$OUT_DIR"
OUT_JSON="${OUT_DIR}/pt1-report.json"

# Run the Python reporter
python3 "$PY_TOOL" \
	--base-dir "$BASE_DIR" \
	--out "$OUT_JSON"

echo "[OK] PT-1 JSON: $OUT_JSON"
