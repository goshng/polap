#!/usr/bin/env bash
#
# polap-bash-parse-himt-mito-report.sh
# Version: 0.3.0
#
# Parse a HiMT mitochondrial HTML report and extract:
#   - gene integrity heatmap -> TSV + gene-set completeness summary
#   - basic mitogenome info table (length, contig count, etc.)
#   - contig table
#   - contig-level 3C metrics (using Python, not AWK)
#
# Usage:
#   polap-bash-parse-himt-mito-report.sh [-p PREFIX] himt_mitochondrial.html
#
# Outputs (PREFIX defaults to input path without .html):
#   PREFIX.gene_integrity.tsv
#   PREFIX.gene_integrity_summary.tsv
#   PREFIX.mito_basic_info.tsv
#   PREFIX.contig_table.tsv
#   PREFIX.contig_3C_metrics.tsv
#

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
	cat <<EOF
Usage: $(basename "$0") [-p PREFIX] himt_mitochondrial.html

Options:
  -p PREFIX   Output prefix (default: input path without .html)
  -h          Show this help message

This script expects the helper Python scripts:
  scripts/himt_parse_mito_html.py
  scripts/contig_3C_metrics.py
relative to the location of this bash script.
EOF
}

prefix=""
while getopts ":p:h" opt; do
	case "${opt}" in
	p)
		prefix="${OPTARG}"
		;;
	h)
		usage
		exit 0
		;;
	\?)
		echo "ERROR: Invalid option -${OPTARG}" >&2
		usage >&2
		exit 1
		;;
	:)
		echo "ERROR: Option -${OPTARG} requires an argument." >&2
		usage >&2
		exit 1
		;;
	esac
done
shift $((OPTIND - 1))

if [ "$#" -ne 1 ]; then
	echo "ERROR: Exactly one HTML input file is required." >&2
	usage >&2
	exit 1
fi

html="$1"

if [ ! -f "$html" ]; then
	echo "ERROR: HTML file not found: $html" >&2
	exit 1
fi

if [ -z "$prefix" ]; then
	# Default prefix: strip trailing .html (if present)
	prefix="${html%.html}"
fi

if ! command -v python3 >/dev/null 2>&1; then
	echo "ERROR: python3 is required but not found in PATH." >&2
	exit 1
fi

parse_py="${script_dir}/scripts/himt_parse_mito_html.py"
contig_metrics_py="${script_dir}/scripts/contig_3C_metrics.py"

if [ ! -s "$parse_py" ]; then
	echo "ERROR: Helper script not found or not executable: $parse_py" >&2
	exit 1
fi

if [ ! -s "$contig_metrics_py" ]; then
	echo "ERROR: Helper script not found or not executable: $contig_metrics_py" >&2
	exit 1
fi

# 1) Extract tables and gene integrity from the HTML
python3 "$parse_py" "$html" "$prefix"

# 2) Compute contig-level 3C metrics using Python (replacement for AWK script)
contig_table="${prefix}.contig_table.tsv"
contig_metrics="${prefix}.contig_3C_metrics.tsv"

if [ ! -f "$contig_table" ]; then
	echo "ERROR: Expected contig table not found: $contig_table" >&2
	exit 1
fi

python3 "$contig_metrics_py" "$contig_table" "$contig_metrics"

cat <<EOF
[polap-bash-parse-himt-mito-report.sh] Done.

Generated files (prefix: $prefix):
  - ${prefix}.gene_integrity.tsv
  - ${prefix}.gene_integrity_summary.tsv
  - ${prefix}.mito_basic_info.tsv
  - ${prefix}.contig_table.tsv
  - ${prefix}.contig_3C_metrics.tsv
EOF
