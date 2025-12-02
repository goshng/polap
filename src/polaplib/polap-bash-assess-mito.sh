#!/usr/bin/env bash
# polap-bash-assess-mito.sh
# Version: v0.1.0
#
# High-level mito assembly assessment pipeline:
#   1) Parse HiMT mitochondrial HTML report -> core TSVs
#   2) Compute 4C metrics -> PREFIX.4c_metrics.tsv
#   3) Generate 4C QC dashboard PDF
#
# Usage:
#   polap-bash-assess-mito.sh -i himt_mitochondrial.html -o PREFIX
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
	cat <<EOF
Usage: $(basename "$0") -i himt_mitochondrial.html -o PREFIX

Options:
  -i FILE   HiMT mitochondrial HTML report
  -o PREFIX Output prefix (e.g. path/to/sample.mito)
  -h        Show this help

Generated files:
  PREFIX.gene_integrity.tsv
  PREFIX.gene_integrity_summary.tsv
  PREFIX.spliced_genes.tsv
  PREFIX.mito_basic_info.tsv
  PREFIX.contig_table.tsv
  PREFIX.4c_metrics.tsv
  PREFIX.mito_qc_4c.pdf
EOF
}

INPUT_HTML=""
PREFIX=""

while getopts ":i:o:h" opt; do
	case "$opt" in
	i) INPUT_HTML="$OPTARG" ;;
	o) PREFIX="$OPTARG" ;;
	h)
		usage
		exit 0
		;;
	\?)
		echo "ERROR: Invalid option -$OPTARG" >&2
		usage >&2
		exit 1
		;;
	:)
		echo "ERROR: Option -$OPTARG requires argument." >&2
		usage >&2
		exit 1
		;;
	esac
done
shift $((OPTIND - 1))

if [[ -z "$INPUT_HTML" || -z "$PREFIX" ]]; then
	echo "ERROR: -i and -o are required." >&2
	usage >&2
	exit 1
fi

if [[ ! -f "$INPUT_HTML" ]]; then
	echo "ERROR: Input HTML not found: $INPUT_HTML" >&2
	exit 1
fi

# paths to helper scripts
PY_PARSE="${SCRIPT_DIR}/scripts/himt_parse_mito_html.py"
PY_4C="${SCRIPT_DIR}/scripts/mt_mito_4c_metrics.py"
R_QC="${SCRIPT_DIR}/scripts/mito_qc_dashboard_4c.R"

command -v python3 >/dev/null 2>&1 || {
	echo "ERROR: python3 not in PATH." >&2
	exit 1
}
command -v Rscript >/dev/null 2>&1 || {
	echo "ERROR: Rscript not in PATH." >&2
	exit 1
}

[[ -s "$PY_PARSE" ]] || {
	echo "ERROR: $PY_PARSE not executable." >&2
	exit 1
}
[[ -s "$PY_4C" ]] || {
	echo "ERROR: $PY_4C not executable." >&2
	exit 1
}
[[ -s "$R_QC" ]] || {
	echo "ERROR: $R_QC not executable." >&2
	exit 1
}

# 1) Parse HTML -> core TSVs
python3 "$PY_PARSE" "$INPUT_HTML" "$PREFIX"

# 2) Compute 4C metrics
python3 "$PY_4C" "$PREFIX"

# 3) QC dashboard
Rscript "$R_QC" \
	"${PREFIX}.4c_metrics.tsv" \
	"${PREFIX}.contig_table.tsv" \
	"${PREFIX}.spliced_genes.tsv" \
	"${PREFIX}.mito_qc_4c.pdf"

echo "[polap-bash-assess-mito] Done."
echo "  4C metrics TSV : ${PREFIX}.4c_metrics.tsv"
echo "  QC dashboard   : ${PREFIX}.mito_qc_4c.pdf"
