#!/usr/bin/env bash
# polap-bash-oatk-ptmt-sheet.sh
# Version : v0.3.0  (2025-12-08)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Build a PT/MT graph sheet PDF for all OATK patterns.
# This script now:
#   1) Runs polap-py-oatk-graph-table.py to generate a pattern/cov PNG CSV
#   2) Runs polap-r-oatk-ptmt-sheet.R to plot one page per pattern

set -euo pipefail

usage() {
	cat <<EOF
Usage: $(basename "$0") [options] [extra R options]

Options:
  --manifest PATH   Input manifest JSON
                    (default: \$POLAP_OATK_MANIFEST or oatk-manifest.json)

  --csv PATH        Output intermediate CSV containing:
                        species,code2,pattern,cov,mito_png,plastid_png
                    (default: \$POLAP_OATK_PTMTSHEET_LIST or oatk-ptmt-sheet-list.csv)

  --out PATH        Output PDF (default: \$POLAP_OATK_PTMTSHEET_PDF or oatk-ptmt-sheet.pdf)
  --pdf PATH        Alias for --out

  -h, --help        Show this help

All extra arguments after options are passed directly to the R script
(polap-r-oatk-ptmt-sheet.R), e.g.:
  --rows-per-page 8 --page-width-in 10 --page-height-in 14
EOF
	exit 1
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
PYTHON="${PYTHON:-python3}"
R_BIN="${R_BIN:-Rscript}"

PY_GRAPH_TABLE="${SCRIPT_DIR}/scripts/polap-py-oatk-graph-table.py"
R_PTMTSHEET="${SCRIPT_DIR}/scripts/polap-r-oatk-ptmt-sheet.R"

MANIFEST="${POLAP_OATK_MANIFEST:-oatk-manifest.json}"
CSV="${POLAP_OATK_PTMTSHEET_LIST:-oatk-ptmt-sheet-list.csv}"
OUT="${POLAP_OATK_PTMTSHEET_PDF:-oatk-ptmt-sheet.pdf}"

R_EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
	case "$1" in
	--manifest)
		MANIFEST="$2"
		shift 2
		;;
	--csv)
		CSV="$2"
		shift 2
		;;
	--out | --pdf)
		OUT="$2"
		shift 2
		;;
	-h | --help)
		usage
		;;
	--)
		shift
		R_EXTRA_ARGS+=("$@")
		break
		;;
	*)
		R_EXTRA_ARGS+=("$1")
		shift
		;;
	esac
done

# ---------------------------------------------------------------------------
# Sanity checks
# ---------------------------------------------------------------------------
if [[ ! -f "$MANIFEST" ]]; then
	echo "[ERR] Manifest JSON not found: $MANIFEST" >&2
	exit 1
fi

mkdir -p "$(dirname "$CSV")" "$(dirname "$OUT")"

# ---------------------------------------------------------------------------
# Step 1 — Build CSV via Python (pattern/cov PNG table)
# ---------------------------------------------------------------------------
echo "[INFO] Generating PT/MT sheet list CSV:"
echo "       $CSV"
"$PYTHON" "$PY_GRAPH_TABLE" \
	--manifest "$MANIFEST" \
	--out-csv "$CSV"

# ---------------------------------------------------------------------------
# Step 2 — Produce multi-page sheet PDF via R
# ---------------------------------------------------------------------------
echo "[INFO] Generating PT/MT sheet PDF:"
echo "       $OUT"

"$R_BIN" "$R_PTMTSHEET" \
	--list "$CSV" \
	--out "$OUT" \
	"${R_EXTRA_ARGS[@]}"

echo "----------------------------------------------------"
echo "Wrote:"
echo "  CSV : $CSV"
echo "  PDF : $OUT"
