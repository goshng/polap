#!/usr/bin/env bash
# File: polap-bash-plant-mt-assembly-trend.sh
# Version: v0.1.0
#
# Build a trend figure of sequencing + assembly methods
# used for plant mitochondrial genomes in NCBI.
#
# INPUT:
#   - Accessions list (one per line), e.g. plant_mt.acc.list
#
# REQUIREMENTS:
#   - Entrez Direct (esearch, efetch)
#   - Python3 + Biopython
#   - Rscript + ggplot2
#
# OUTPUT (under OUTDIR):
#   plant_mt.gbff
#   plant_mt.assembly_meta.tsv
#   plant_mt.assembly_trend_counts.tsv
#   plant_mt.assembly_trend.pdf
#
set -Eeuo pipefail

VERSION="v0.1.0"

usage() {
	cat <<EOF
polap-bash-plant-mt-assembly-trend.sh ${VERSION}
Build trend figure for sequencing and assembly methods used
in published plant mitochondrial genomes (NCBI).

USAGE:
  $(basename "$0") -i plant_mt.acc.list -o outdir

OPTIONS:
  -i FILE      Accession list (one accession per line)
  -o DIR       Output directory
  --batch N    Fetch GenBank in batches of N accessions (default: 100)
  -h, --help   Show this help message

NOTES:
  - Uses efetch -db nucleotide -format gbwithparts to download GenBank.
  - Parses Assembly Method and Sequencing Technology from COMMENT /
    Assembly-Data structured comments via Python.
  - R script then produces a PDF trend figure.
EOF
}

IN_ACC=""
OUTDIR=""
BATCH=100

while [[ $# -gt 0 ]]; do
	case "$1" in
	-i)
		IN_ACC="$2"
		shift 2
		;;
	-o)
		OUTDIR="$2"
		shift 2
		;;
	--batch)
		BATCH="$2"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "ERROR: Unknown argument: $1" >&2
		usage
		exit 1
		;;
	esac
done

if [[ -z "$IN_ACC" || -z "$OUTDIR" ]]; then
	echo "ERROR: -i and -o are required" >&2
	usage
	exit 1
fi

mkdir -p "$OUTDIR"

need_cmd() {
	command -v "$1" >/dev/null 2>&1 || {
		echo "ERROR: Required command not found in PATH: $1" >&2
		exit 1
	}
}

need_cmd esearch
need_cmd efetch
need_cmd python3
need_cmd Rscript

GBFF="${OUTDIR}/plant_mt.gbff"
META_TSV="${OUTDIR}/plant_mt.assembly_meta.tsv"
COUNTS_TSV="${OUTDIR}/plant_mt.assembly_trend_counts.tsv"
PLOT_PDF="${OUTDIR}/plant_mt.assembly_trend.pdf"

# 1) Download GenBank (batched)
echo "[1/3] Downloading GenBank records to: $GBFF" >&2
: >"$GBFF"

batch_fetch() {
	local batch_file="$1"
	local ids
	ids=$(paste -sd, "$batch_file")
	if [[ -z "$ids" ]]; then
		return
	fi
	efetch -db nucleotide -id "$ids" -format gbwithparts >>"$GBFF"
}

tmp_batch="$(mktemp)"
trap 'rm -f "$tmp_batch"' EXIT

i=0
: >"$tmp_batch"
while read -r acc; do
	[[ -z "$acc" ]] && continue
	echo "$acc" >>"$tmp_batch"
	((i++))
	if ((i >= BATCH)); then
		batch_fetch "$tmp_batch"
		: >"$tmp_batch"
		i=0
	fi
done <"$IN_ACC"

# fetch remaining
if [[ -s "$tmp_batch" ]]; then
	batch_fetch "$tmp_batch"
fi

echo "  GenBank download complete." >&2

# 2) Parse assembly + sequencing metadata with Python
echo "[2/3] Parsing assembly and sequencing metadata ..." >&2

python3 scripts/parse_mt_assembly_from_genbank.py \
	"$GBFF" \
	"$META_TSV"

echo "  Parsed metadata written: $META_TSV" >&2

# 3) Plot trend with R
echo "[3/3] Creating trend figure with R ..." >&2

Rscript scripts/plot_mt_assembly_trend.R \
	--meta "$META_TSV" \
	--out-pdf "$PLOT_PDF" \
	--out-tsv "$COUNTS_TSV"

echo "  Trend counts written: $COUNTS_TSV" >&2
echo "  Trend figure written: $PLOT_PDF" >&2
echo "Done." >&2
