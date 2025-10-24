#!/usr/bin/env bash
# polap-bash-report-readassemble-pt-json2html.sh
# Version: v0.3.0
# SPDX-License-Identifier: GPL-3.0-or-later
#
# Convert pt1-report.json -> HTML using a separate Python script.
#
# Usage:
#   polap-bash-report-readassemble-pt-json2html.sh --json /path/to/pt1-report.json --html /path/to/pt1-report.html
#
set -euo pipefail
IFS=$'\n\t'

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --json PT1_REPORT.JSON --html OUT.HTML
EOF
}

JSON_IN=""
HTML_OUT=""

while (($#)); do
	case "$1" in
	--json)
		JSON_IN="${2:?}"
		shift 2
		;;
	--html)
		HTML_OUT="${2:?}"
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

[[ -n "$JSON_IN" && -f "$JSON_IN" ]] || {
	echo "[ERR] --json missing or not a file: $JSON_IN" >&2
	exit 2
}
[[ -n "$HTML_OUT" ]] || {
	echo "[ERR] --html missing" >&2
	exit 2
}
mkdir -p "$(dirname "$HTML_OUT")"

# Locate sibling python tool
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POLAPLIB_DIR="${POLAPLIB_DIR:-$SCRIPT_DIR}"
PY_TOOL="${POLAPLIB_DIR}/scripts/report_readassemble_pt_json2html.py"

[[ -f "$PY_TOOL" ]] || {
	echo "[ERR] Missing Python helper: $PY_TOOL" >&2
	exit 3
}

python3 "$PY_TOOL" --json "$JSON_IN" --html "$HTML_OUT"

echo "[OK] HTML: $HTML_OUT"
