#!/usr/bin/env bash
# File: polap-bash-crawl-plant-mtgenomes.sh
# Version: v0.1.0
#
# Crawl NCBI for plant mitochondrial complete genomes (GenBank)
# and associated SRA runs (RunInfo CSV).
#
# REQUIREMENTS:
#   - Entrez Direct (esearch, efetch, xtract)
#     https://www.ncbi.nlm.nih.gov/books/NBK179288/  (Entrez Direct)
#   - (optional) SRA Toolkit for downloading reads (prefetch, fasterq-dump)
#
# OUTPUT (under OUTDIR):
#   plant_mt.query.txt        : the exact Entrez query used
#   plant_mt.acc.list         : accession list (one per line)
#   plant_mt.metadata.tsv     : basic metadata per accession
#   plant_mt.sra_runinfo.csv  : aggregated SRA RunInfo (if --with-sra)
#
# EXAMPLE:
#   polap-bash-crawl-plant-mtgenomes.sh -o man/analysis/plant_mt_crawl --with-sra
#
set -Eeuo pipefail

VERSION="v0.1.0"

usage() {
	cat <<EOF
polap-bash-crawl-plant-mtgenomes.sh ${VERSION}
Crawl NCBI for plant mitochondrial complete genomes (GenBank + SRA metadata).

USAGE:
  $(basename "$0") -o <outdir> [--with-sra] [--max N] [--query "custom Entrez query"]

OPTIONS:
  -o DIR         Output directory (required)
  --with-sra     Also query SRA to build plant_mt.sra_runinfo.csv
  --max N        Limit number of nucleotide records (for testing)
  --query Q      Override default Entrez query (see below)
  -h, --help     Show this help message and exit

DEFAULT QUERY (db=nucleotide):
  (mitochondrion[filter]) AND plants[organism] AND (\"complete genome\"[Title] OR \"complete sequence\"[Title])

NCBI QUERY:
  mitochondrion[filter] AND "Embryophyta"[Organism] AND ("complete genome"[Title] OR "complete sequence"[Title])

NOTES:
  - Requires Entrez Direct (esearch/efetch/xtract).
  - Uses efetch -format runinfo to get SRA RunInfo CSV for each accession
    as described in NCBI documentation.
EOF
}

# ---------- Args parsing ----------
OUTDIR=""
WITH_SRA=0
MAX_REC=""
CUSTOM_QUERY=""
SKIP_NT=0

while [[ $# -gt 0 ]]; do
	case "$1" in
	-o)
		OUTDIR="$2"
		shift 2
		;;
	--with-sra)
		WITH_SRA=1
		shift
		;;
	--max)
		MAX_REC="$2"
		shift 2
		;;
	--query)
		CUSTOM_QUERY="$2"
		shift 2
		;;
	--skip-nt)
		SKIP_NT=1
		shift
		;; # NEW
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

if [[ -z "$OUTDIR" ]]; then
	echo "ERROR: -o <outdir> is required" >&2
	usage
	exit 1
fi

mkdir -p "$OUTDIR"

# ---------- Check dependencies ----------
need_cmd() {
	command -v "$1" >/dev/null 2>&1 || {
		echo "ERROR: Required command not found in PATH: $1" >&2
		exit 1
	}
}

need_cmd esearch
need_cmd efetch
need_cmd xtract

# ---------- Build query ----------
DEFAULT_QUERY='(mitochondrion[filter]) AND plants[organism] AND ("complete genome"[Title] OR "complete sequence"[Title])'

QUERY="${CUSTOM_QUERY:-$DEFAULT_QUERY}"

echo "Query: $QUERY"
echo "$QUERY" >"${OUTDIR}/plant_mt.query.txt"

# ---------- Step 1: get accession list ----------
ACC_LIST="${OUTDIR}/plant_mt.acc.list"
META_TSV="${OUTDIR}/plant_mt.metadata.tsv"
RUNINFO_CSV="${OUTDIR}/plant_mt.sra_runinfo.csv"

if ((SKIP_NT == 0)); then
	# ---------- Step 1: get accession list ----------
	echo "[1/3] Fetching accession list from NCBI nucleotide ..." >&2

	if [[ -n "$MAX_REC" ]]; then
		esearch -db nucleotide -query "$QUERY" -retmax "$MAX_REC" |
			efetch -format acc \
				>"$ACC_LIST"
	else
		esearch -db nucleotide -query "$QUERY" |
			efetch -format acc \
				>"$ACC_LIST"
	fi

	N_ACC=$(wc -l <"$ACC_LIST" || echo 0)
	echo "  Found $N_ACC nucleotide accessions" >&2

	if ((N_ACC == 0)); then
		echo "WARNING: No accessions found. Check your query." >&2
		exit 0
	fi

	# ---------- Step 2: fetch basic metadata (tsv) ----------
	echo "[2/3] Fetching nucleotide metadata (esummary + xtract) ..." >&2

	if [[ -n "$MAX_REC" ]]; then
		esearch -db nucleotide -query "$QUERY" -retmax "$MAX_REC" |
			esummary |
			xtract -pattern DocumentSummary \
				-element AccessionVersion Organism Title Length CreateDate BioProject \
				>"$META_TSV"
	else
		esearch -db nucleotide -query "$QUERY" |
			esummary |
			xtract -pattern DocumentSummary \
				-element AccessionVersion Organism Title Length CreateDate BioProject \
				>"$META_TSV"
	fi

	tmp_meta="${META_TSV}.tmp"
	{
		echo -e "AccessionVersion\tOrganism\tTitle\tLength\tCreateDate\tBioProject"
		cat "$META_TSV"
	} >"$tmp_meta"
	mv "$tmp_meta" "$META_TSV"

	echo "  Metadata written: $META_TSV" >&2

else
	echo "[1/3] SKIP_NT=1: skipping accession+metadata fetch" >&2
	# Sanity checks
	if [[ ! -s "$ACC_LIST" ]]; then
		echo "ERROR: --skip-nt was set but $ACC_LIST does not exist or is empty" >&2
		exit 1
	fi
	if [[ ! -s "$META_TSV" ]]; then
		echo "WARNING: --skip-nt was set but $META_TSV does not exist or is empty" >&2
	fi
fi

# ---------- Step 3: optional SRA RunInfo aggregation ----------
if ((WITH_SRA == 1)); then
	echo "[3/3] Fetching SRA RunInfo for each accession via nuccore->SRA links ..." >&2
	need_cmd elink

	: >"$RUNINFO_CSV"
	header_done=0

	while read -r ACC; do
		[[ -z "$ACC" ]] && continue
		echo "  [SRA] $ACC" >&2

		RUNINFO_TMP="$(mktemp)"
		set +e
		# Use nuccore->sra link instead of searching SRA by string
		elink -db nucleotide -id "$ACC" -target sra |
			efetch -format runinfo >"$RUNINFO_TMP" 2>/dev/null
		# elink -db nuccore -target sra -id "$ACC" |
		# 	efetch -db sra -format runinfo >"$RUNINFO_TMP" 2>/dev/null
		rc=$?
		set -e

		if ((rc != 0)) || [[ ! -s "$RUNINFO_TMP" ]]; then
			rm -f "$RUNINFO_TMP"
			continue
		fi

		if ((header_done == 0)); then
			cat "$RUNINFO_TMP" >>"$RUNINFO_CSV"
			header_done=1
		else
			# append but skip header line
			awk 'NR>1' "$RUNINFO_TMP" >>"$RUNINFO_CSV"
		fi
		rm -f "$RUNINFO_TMP"
	done <"$ACC_LIST"

	if [[ -s "$RUNINFO_CSV" ]]; then
		echo "  SRA RunInfo written: $RUNINFO_CSV" >&2
	else
		echo "  WARNING: No SRA RunInfo links found for any accession." >&2
	fi
else
	echo "[3/3] Skipping SRA RunInfo (no --with-sra)" >&2
fi

echo "Done. Outputs in: $OUTDIR" >&2
