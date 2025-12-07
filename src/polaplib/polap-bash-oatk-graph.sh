#!/usr/bin/env bash
# polap-bash-oatk-graph.sh
# Version : v0.1.0  (2025-12-07)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Build long-format metrics CSV + QC boxplot PDF for Oatk (and other pipelines)
# from an auto-manifest JSON.
#
# Defaults from env:
#   POLAP_OATK_MANIFEST
#   POLAP_OATK_GRAPH_CSV
#   POLAP_OATK_GRAPH_PDF
#   POLAP_OATK_GRAPH_TITLE
#
# CLI:
#   --manifest PATH   (default: \$POLAP_OATK_MANIFEST or oatk-manifest.json)
#   --csv PATH        (default: \$POLAP_OATK_GRAPH_CSV or oatk-graph-metrics.csv)
#   --pdf PATH        (default: \$POLAP_OATK_GRAPH_PDF or oatk-graph-metrics.pdf)
#   --title STRING    (default: env or script default title)
#
# Any other arguments are forwarded to the R script, e.g.:
#   --ncol-per-row 2 --page-width-in 8 --page-height-in 6
#

set -euo pipefail

usage() {
	cat <<EOF
Usage: $(basename "$0") [options] [extra R options]

Options:
  --manifest PATH   Input manifest JSON (default: \$POLAP_OATK_MANIFEST or oatk-manifest.json)
  --csv PATH        Output metrics CSV (default: \$POLAP_OATK_GRAPH_CSV or oatk-graph-metrics.csv)
  --pdf PATH        Output PDF (default: \$POLAP_OATK_GRAPH_PDF or oatk-graph-metrics.pdf)
  --title STRING    Plot title (default: \$POLAP_OATK_GRAPH_TITLE or built-in default)
  -h, --help        Show this help

Any remaining arguments are passed through to polap-r-oatk-graph.R.
EOF
	exit 1
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
PYTHON="${PYTHON:-python3}"
R_BIN="${R_BIN:-Rscript}"

GRAPH_TABLE_PY="${SCRIPT_DIR}/scripts/polap-py-oatk-graph-table.py"
GRAPH_R="${SCRIPT_DIR}/scripts/polap-r-oatk-graph.R"

MANIFEST="${POLAP_OATK_MANIFEST:-oatk-manifest.json}"
CSV="${POLAP_OATK_GRAPH_CSV:-oatk-graph-metrics.csv}"
PDF="${POLAP_OATK_GRAPH_PDF:-oatk-graph-metrics.pdf}"
TITLE="${POLAP_OATK_GRAPH_TITLE:-Mitogenome assembly QC metrics (PMAT, TIPPo, HiMT, Oatk)}"

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
	--pdf | --out)
		PDF="$2"
		shift 2
		;;
	--title)
		TITLE="$2"
		shift 2
		;;
	-h | --help) usage ;;
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

if [[ ! -f "$MANIFEST" ]]; then
	echo "ERROR: manifest JSON not found: $MANIFEST" >&2
	exit 1
fi

mkdir -p "$(dirname "$CSV")" "$(dirname "$PDF")"

# 1) Manifest -> long-format CSV
"$PYTHON" "$GRAPH_TABLE_PY" \
	--manifest "$MANIFEST" \
	--out-csv "$CSV"

# 2) CSV -> boxplot PDF
"$R_BIN" "$GRAPH_R" \
	--data "$CSV" \
	--out "$PDF" \
	--title "$TITLE" \
	"${R_EXTRA_ARGS[@]}"

echo "Wrote:"
echo "  CSV : $CSV"
echo "  PDF : $PDF"
