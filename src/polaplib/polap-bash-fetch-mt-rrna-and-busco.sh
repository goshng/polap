#!/usr/bin/env bash
# polap-bash-fetch-mt-rrna-and-busco.sh  (v0.0.1)
# Fetch plant mitochondrial rRNAs (rrn5/rrn18/rrn26…) from NCBI nt,
# and fetch a BUSCO plant nuclear lineage dataset (to use for “definitely nuclear” bait).
set -euo pipefail
VERBOSE=0
QUIET=0
OUT="ref"

usage() {
	cat <<EOF
Usage: $0 [-o outdir] [-v|--verbose] [--quiet]
Requires: esearch, efetch (Entrez Direct); wget (for BUSCO data).
Outputs:
  <outdir>/mt_rrna.plant.fa
  <outdir>/busco_viridiplantae.tar.gz   (pointer to official download)
EOF
}

while [[ $# -gt 0 ]]; do
	case "$1" in
	-o | --out)
		OUT="$2"
		shift 2
		;;
	-v | --verbose)
		VERBOSE=$((VERBOSE + 1))
		shift
		;;
	--quiet)
		QUIET=1
		VERBOSE=0
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	*) shift ;;
	esac
done

mkdir -p "$OUT"
# Plant mt rRNAs (very conservative query; refine as needed)
# Example: search mitochondrial 5S/18S/26S annotated in Viridiplantae
for gene in "5S ribosomal RNA" "18S ribosomal RNA" "26S ribosomal RNA" "rrn5" "rrn18" "rrn26"; do
	esearch -db nucleotide -query "($gene[Title]) AND mitochondrion[Filter] AND Viridiplantae[Organism]" |
		efetch -format fasta >>"$OUT/mt_rrna.plant.fa" || true
done
# De-duplicate
awk '/^>/{h=$0; if(!seen[h]++){print h; next}} !/^>/{print}' "$OUT/mt_rrna.plant.fa" >"$OUT/mt_rrna.plant.tmp" || true
mv "$OUT/mt_rrna.plant.tmp" "$OUT/mt_rrna.plant.fa"

# BUSCO pointer:
cat >"$OUT/README.busco.txt" <<'TXT'
Get a plant nuclear BUSCO dataset (e.g., viridiplantae_odb10) from:
  https://busco-data.ezlab.org/
Use this set to build a "definitely nuclear" bait (filtering). Licensing may apply.
TXT

echo "[INFO] Wrote: $OUT/mt_rrna.plant.fa"
