#!/usr/bin/env bash
# File: polap-bash-mt-meta-analysis.sh
# Version: v0.1.0
#
# Orchestrate plant mitochondrial genome meta-analysis:
#  1. Fetch GBFF for all accessions (if missing)
#  2. Parse GBFF to assembly_meta.tsv (Python)
#  3. Collect SRA RunInfo via BioProject -> SRA (EDirect)
#  4. Generate overview plots & summary TSVs (R)
#
# USAGE:
#   polap-bash-mt-meta-analysis.sh -i plant_mt.acc.list -o outdir
#
# REQUIREMENTS:
#   - EDirect installed under $HOME/edirect (esearch, efetch)
#   - Python 3 + Biopython
#   - R + ggplot2 + dplyr + readr
#
set -Eeuo pipefail

VERSION="v0.1.0"
POLAPLIB_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"

usage() {
	cat <<EOF
polap-bash-mt-meta-analysis.sh ${VERSION}
Run plant mitochondrial genome meta-analysis (methods + data availability).

USAGE:
  $(basename "$0") -i <accession_list> -o <outdir> [--batch N] [--force-gbff]

OPTIONS:
  -i FILE      Accession list (one accession per line)
  -o DIR       Output directory
  --batch N    Fetch GBFF in batches of N accessions (default: 100)
  --force-gbff Re-download GBFF even if it already exists
  -h, --help   Show this help message and exit

STEPS:
  1) Fetch GenBank/GBFF for all mitogenome accessions
  2) Parse GBFF -> assembly_meta.tsv (methods, platform, BioProject, SRA)
  3) BioProject -> SRA RunInfo (project-level)
  4) Plot trends & data availability (PDFs + TSVs)
EOF
}

IN_ACC=""
OUTDIR=""
BATCH=100
FORCE_GBFF=0

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
	--force-gbff)
		FORCE_GBFF=1
		shift
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

[[ -z "$IN_ACC" || -z "$OUTDIR" ]] && {
	echo "ERROR: -i and -o are required" >&2
	usage
	exit 1
}

mkdir -p "$OUTDIR"

# ----- EDirect paths -----
: "${EDIRECT_DIR:=$HOME/edirect}"
ESEARCH="esearch"
EFETCH="efetch"

# if [[ ! -x "$ESEARCH" || ! -x "$EFETCH" ]]; then
# 	echo "ERROR: EDirect esearch/efetch not found under ${EDIRECT_DIR}" >&2
# 	echo "       Install Entrez Direct or set EDIRECT_DIR appropriately." >&2
# 	exit 1
# fi

# ----- Check Python/R -----
command -v python3 >/dev/null 2>&1 || {
	echo "ERROR: python3 not in PATH" >&2
	exit 1
}
command -v Rscript >/dev/null 2>&1 || {
	echo "ERROR: Rscript not in PATH" >&2
	exit 1
}

GBFF="${OUTDIR}/plant_mt.gbff"
META_TSV="${OUTDIR}/plant_mt.assembly_meta.tsv"
SRA_RUNINFO_CSV="${OUTDIR}/plant_mt.sra_runinfo_by_bioproject.csv"

echo "[meta] IN_ACC = $IN_ACC"
echo "[meta] OUTDIR = $OUTDIR"
echo "[meta] GBFF   = $GBFF"

# ============================================================================
# STEP 1: Fetch GBFF
# ============================================================================
if [[ -s "$GBFF" && $FORCE_GBFF -eq 0 ]]; then
	echo "[1/4] GBFF already exists (use --force-gbff to re-download): $GBFF"
else
	echo "[1/4] Fetching GenBank/GBFF for all accessions into: $GBFF"
	>"$GBFF"

	tmp_batch="$(mktemp)"
	trap 'rm -f "$tmp_batch"' EXIT

	i=0
	while read -r acc; do
		[[ -z "$acc" ]] && continue
		echo "$acc" >>"$tmp_batch"
		i=$((i + 1))
		if ((i >= BATCH)); then
			ids=$(paste -sd, "$tmp_batch")
			echo "  [GBFF] fetching batch of ${i} accessions..." >&2
			"$EFETCH" -db nucleotide -id "$ids" -format gbwithparts >>"$GBFF"
			: >"$tmp_batch"
			i=0
		fi
	done <"$IN_ACC"

	# flush last partial batch
	if [[ -s "$tmp_batch" ]]; then
		ids=$(paste -sd, "$tmp_batch")
		echo "  [GBFF] fetching final batch of ${i} accessions..." >&2
		"$EFETCH" -db nucleotide -id "$ids" -format gbwithparts >>"$GBFF"
	fi

	rm -f "$tmp_batch"
	echo "  [GBFF] written: $GBFF"
fi

# ============================================================================
# STEP 2: Parse GBFF -> assembly_meta.tsv (Python)
# ============================================================================
echo "[2/4] Parsing GBFF to assembly_meta.tsv ..."
python3 "${POLAPLIB_DIR}/scripts/parse_mt_assembly_from_genbank.py" \
	"$GBFF" \
	"$META_TSV"

echo "  [META] written: $META_TSV"

# ============================================================================
# STEP 3: BioProject -> SRA RunInfo via EDirect
# ============================================================================
echo "[3/4] Collecting SRA RunInfo via BioProject -> SRA links ..."

# Extract unique BioProjects (column 8 in META_TSV)
BIOPROJ_LIST="${OUTDIR}/plant_mt.bioproject.list"
awk -F'\t' 'NR>1 && $8 != "" {print $8}' "$META_TSV" | sort -u >"$BIOPROJ_LIST"

N_BP=$(wc -l <"$BIOPROJ_LIST" || echo 0)
echo "  [SRA] Found ${N_BP} unique BioProject IDs in assembly_meta.tsv"

>"$SRA_RUNINFO_CSV"
header_done=0

while read -r PRJ; do
	[[ -z "$PRJ" ]] && continue
	echo "  [SRA] BioProject: $PRJ" >&2

	tmp_csv="$(mktemp)"
	set +e
	"$ESEARCH" -db sra -query "$PRJ" | "$EFETCH" -format runinfo >"$tmp_csv" 2>/dev/null
	rc=$?
	set -e

	if ((rc != 0)) || [[ ! -s "$tmp_csv" ]]; then
		rm -f "$tmp_csv"
		continue
	fi

	if ((header_done == 0)); then
		cat "$tmp_csv" >>"$SRA_RUNINFO_CSV"
		header_done=1
	else
		# append without header
		awk 'NR>1' "$tmp_csv" >>"$SRA_RUNINFO_CSV"
	fi
	rm -f "$tmp_csv"
done <"$BIOPROJ_LIST"

if [[ -s "$SRA_RUNINFO_CSV" ]]; then
	echo "  [SRA] SRA RunInfo written: $SRA_RUNINFO_CSV"
else
	echo "  [SRA] WARNING: No SRA RunInfo found for any BioProject."
fi

# ============================================================================
# STEP 4: Plot overview trends & summaries (R)
# ============================================================================
echo "[4/4] Generating overview plots and summary tables (R) ..."
Rscript --vanilla "${POLAPLIB_DIR}/scripts/plot_mt_meta_overview.R" \
	--meta "$META_TSV" \
	--sra "$SRA_RUNINFO_CSV" \
	--outdir "$OUTDIR"

echo "[meta] All done. Outputs in: $OUTDIR"
