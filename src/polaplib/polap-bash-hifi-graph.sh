#!/usr/bin/env bash
# polap-bash-hifi-graph.sh
# Version : v1.1.1  (2025-12-02)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Build boxplot PDF of mitogenome assembly QC metrics across 4 pipelines:
#   pmat (pmat2), tippo, himt, oatk (c30 config if available)
#
# Metrics (columns in CSV for R):
#   mem_gb, time_hours, disk_gb,
#   geneset_completeness_prop,
#   num_contigs, N50, fragmentation_index,
#   contig_length_cv, max_contig_prop, total_length
#
set -euo pipefail
IFS=$'\n\t'

_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export _POLAPLIB_DIR

# -----------------------------------------------------------------------------#
# Defaults and CLI args
# -----------------------------------------------------------------------------#
BASE_DIR=""
MANIFEST=""
OUTPDF=""
TITLE="Mitogenome assembly QC metrics (PMAT, TIPPo, HiMT, Oatk)"
PW=11 PH=8.5 MG=0.5

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --manifest v5-0-auto-manifest.json --out md/hifi-mt-graph.pdf \\
                   [--base-dir <dir>] \\
                   [--page-width-in IN] [--page-height-in IN] [--margin-in IN] \\
                   [--title "Mitogenome assembly QC metrics (PMAT, TIPPo, HiMT, Oatk)"]

Required:
  --manifest   v5-0-auto-manifest.json for this analysis
  --out        Output PDF path

This script creates a CSV with columns:
  species, pipeline,
  mem_gb, time_hours, disk_gb,
  geneset_completeness_prop,
  num_contigs, N50, fragmentation_index,
  contig_length_cv, max_contig_prop, total_length
which is passed to polap-r-hifi-graph.R.
EOF
}

while (($#)); do
	case "$1" in
	--base-dir)
		BASE_DIR="${2:?}"
		shift 2
		;;
	--manifest)
		MANIFEST="${2:?}"
		shift 2
		;;
	--out)
		OUTPDF="${2:?}"
		shift 2
		;;
	--page-width-in)
		PW="${2:?}"
		shift 2
		;;
	--page-height-in)
		PH="${2:?}"
		shift 2
		;;
	--margin-in)
		MG="${2:?}"
		shift 2
		;;
	--title)
		TITLE="${2:?}"
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

[[ -s "$MANIFEST" && -n "$OUTPDF" ]] || {
	echo "[ERR] need --manifest and --out" >&2
	exit 2
}

mkdir -p "$(dirname "$OUTPDF")"

# -----------------------------------------------------------------------------#
# Determine root (for completeness – not actually needed for metrics)
# -----------------------------------------------------------------------------#
if [[ -n "$BASE_DIR" ]]; then
	SR="$(cd "$(dirname "$BASE_DIR")" && pwd)"
else
	SR="$(cd "$(dirname "$MANIFEST")/.." && pwd)"
fi

# -----------------------------------------------------------------------------#
# Temp files
# -----------------------------------------------------------------------------#
CSV="$(mktemp -t hifi-mt-metrics.XXXXXX.csv)"
trap 'rm -f "$CSV"' EXIT

# -----------------------------------------------------------------------------#
# Build long-format metrics table from manifest
#   One row per (species, pipeline)
#   Pipelines: pmat (pmat2), tippo, himt, oatk (c30_* metrics)
# -----------------------------------------------------------------------------#
if ! command -v jq >/dev/null 2>&1; then
	echo "[ERR] jq required" >&2
	exit 2
fi

{
	echo "species,pipeline,mem_gb,time_hours,disk_gb,geneset_completeness_prop,num_contigs,N50,fragmentation_index,contig_length_cv,max_contig_prop,total_length"
	jq -r '
	  .items[]
	  | [
	      # PMAT (pmat2)
	      [
	        .species, "pmat",
	        (.pmat2.mem_gb // ""),
	        (.pmat2.time_hours // ""),
	        (.pmat2.disk_gb // ""),
	        (.pmat2.geneset_completeness_prop // ""),
	        (.pmat2.num_contigs // ""),
	        (.pmat2.N50 // ""),
	        (.pmat2.fragmentation_index // ""),
	        (.pmat2.contig_length_cv // ""),
	        (.pmat2.max_contig_prop // ""),
	        (.pmat2.total_length // "")
	      ],
	      # TIPPo
	      [
	        .species, "tippo",
	        (.tippo.mem_gb // ""),
	        (.tippo.time_hours // ""),
	        (.tippo.disk_gb // ""),
	        (.tippo.geneset_completeness_prop // ""),
	        (.tippo.num_contigs // ""),
	        (.tippo.N50 // ""),
	        (.tippo.fragmentation_index // ""),
	        (.tippo.contig_length_cv // ""),
	        (.tippo.max_contig_prop // ""),
	        (.tippo.total_length // "")
	      ],
	      # HiMT
	      [
	        .species, "himt",
	        (.himt.mem_gb // ""),
	        (.himt.time_hours // ""),
	        (.himt.disk_gb // ""),
	        (.himt.geneset_completeness_prop // ""),
	        (.himt.num_contigs // ""),
	        (.himt.N50 // ""),
	        (.himt.fragmentation_index // ""),
	        (.himt.contig_length_cv // ""),
	        (.himt.max_contig_prop // ""),
	        (.himt.total_length // "")
	      ],
	      # Oatk: prefer c30_* metrics
	      [
	        .species, "oatk",
	        (.oatk.c30_mem_gb // ""),
	        (.oatk.c30_time_hours // ""),
	        (.oatk.c30_disk_gb // ""),
	        (.oatk.c30_geneset_completeness_prop // ""),
	        (.oatk.c30_num_contigs // ""),
	        (.oatk.c30_N50 // ""),
	        (.oatk.c30_fragmentation_index // ""),
	        (.oatk.c30_contig_length_cv // ""),
	        (.oatk.c30_max_contig_prop // ""),
	        (.oatk.c30_total_length // "")
	      ]
	    ]
	  | .[]
	  | @csv
	' "$MANIFEST"
} >"$CSV"

# -----------------------------------------------------------------------------#
# Render with R
# -----------------------------------------------------------------------------#
Rscript "${_POLAPLIB_DIR}/scripts/polap-r-hifi-graph.R" \
	--data "$CSV" \
	--out "$OUTPDF" \
	--title "$TITLE" \
	--ncol-per-row 3 \
	--page-width-in "$PW" --page-height-in "$PH" --margin-in "$MG"

echo "[OK] Wrote: $OUTPDF"
