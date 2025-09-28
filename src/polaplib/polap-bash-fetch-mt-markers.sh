#!/usr/bin/env bash
# polap-bash-fetch-mt-markers.sh  v0.0.1
# Fetch plant mitochondrial rRNAs from NCBI, and (optionally) BUSCO plant set.
# Requires: NCBI edirect (esearch/efetch/xtract), curl/wget, gzip
# Usage:
#   polap-bash-fetch-mt-markers.sh -o outdir [--busco] [-v|-q]

set -euo pipefail

ver="v0.0.1"
outdir=""
want_busco=0
verbose=0
quiet=0

log1() { [[ $quiet -eq 0 && $verbose -gt 0 ]] && echo "[INFO] $*" >&2; }
log0() { echo "[ERROR] $*" >&2; }

usage() {
	cat <<EOF
$0 $ver
Fetch plant mitochondrial rRNAs (nuccore FASTA) using NCBI edirect.
Optionally point to BUSCO embryophyta set download.

Options:
  -o, --outdir DIR    output directory (required)
      --busco         also print BUSCO Embryophyta download hints
  -v, --verbose       verbose logs
      --quiet         suppress info logs
  -h, --help          show help
EOF
	exit 0
}

# parse CLI
while [[ $# -gt 0 ]]; do
	case "$1" in
	-o | --outdir)
		outdir="$2"
		shift 2
		;;
	--busco)
		want_busco=1
		shift
		;;
	-v | --verbose)
		verbose=$((verbose + 1))
		shift
		;;
	--quiet)
		quiet=1
		verbose=0
		shift
		;;
	-h | --help) usage ;;
	--version)
		echo "$ver"
		exit 0
		;;
	*)
		log0 "Unknown option: $1"
		usage
		;;
	esac
done

[[ -z "$outdir" ]] && {
	log0 "Missing --outdir"
	exit 2
}
mkdir -p "$outdir"

# sanity: edirect available?
if ! command -v esearch >/dev/null 2>&1; then
	log0 "NCBI edirect (esearch/efetch) not found in PATH."
	echo "Install: https://www.ncbi.nlm.nih.gov/books/NBK179288/ (or conda install -c bioconda entrez-direct)" >&2
	exit 2
fi
if ! command -v cd-hit-est >/dev/null 2>&1; then
	log0 "ch-hit not found."
	log0 "https://github.com/weizhongli/cdhit"
	exit 2
fi

# Plant mitochondrial rRNAs (Embryophyta)
# Query: rRNA features from mitochondrial sequences in Embryophyta
# NCBI nuccore filter tokens:
#   mitochondrion[filter]  ;  rRNA[filter]  ;  Embryophyta[Organism]
log1 "Querying NCBI for Embryophyta mitochondrial rRNAs..."
esearch -db nuccore -query '"Embryophyta"[Organism] AND "18S ribosomal RNA"[Title] OR "26S ribosomal RNA"[Title] AND (mitochondrion[filter] AND ("100"[SLEN] : "3000"[SLEN]))' |
	efetch -format fasta >"$outdir/plant_mt_rRNA.fna"

cd-hit-est -i "$outdir/plant_mt_rRNA.fna" "$outdir/polap-mt.rrna.1.fna"

# Also a stricter query only for the main 5S/18S/26S labels (kept as a subset):
log1 "Querying targeted labels (optional subset): rrn5|5S, rrn18|18S, rrn26|26S..."

# esearch -db nuccore -query '(mitochondrion[filter]) AND (Embryophyta[Organism]) AND (("5S ribosomal RNA"[Title]) OR ("18S ribosomal RNA"[Title]) OR ("26S ribosomal RNA"[Title]) OR rrn5[All Fields] OR rrn18[All Fields] OR rrn26[All Fields])' |
# 	efetch -format fasta >"$outdir/plant_mt_rRNA_core.fna"
#
echo "[OK] Wrote:"
echo "  - $outdir/plant_mt_rRNA.fna"
# echo "  - $outdir/plant_mt_rRNA_core.fna"

if [[ $want_busco -eq 1 ]]; then
	cat <<'NOTE'
[BUSCO Embryophyta (odb10)]
- BUSCO datasets are distributed from the BUSCO/OrthoDB mirrors (license applies).
- If you have a license/permission, download the lineage dataset, e.g.:
    # Example (adjust to current BUSCO/ODB URLs)
    # wget -O embryophyta_odb10.tar.gz "https://busco-data.ezlab.org/v5/data/lineages/embryophyta_odb10.2023-02-23.tar.gz"
    # tar -xzf embryophyta_odb10.tar.gz -C busco_lineages/
- You can then build a “nuclear only” k-mer screen ref using the provided HMM/FASTA markers
  or fetch EnsemblPlants CDS/proteins for your clade as a broader nuclear screen.
NOTE
fi
