#!/usr/bin/env bash
# polap-bash-download-anchors.sh
# Version: v0.1.0
# Description:
#   Download canonical plastid anchor genes (psbA, trnH-GUG, rbcL, matK, ycf1, ndhF)
#   from the Arabidopsis thaliana chloroplast genome (RefSeq NC_000932.1)
#   using NCBI Entrez Direct (efetch) or datasets CLI.
#
# Requirements:
#   conda install -c bioconda entrez-direct seqkit
#   (optional) conda install -c conda-forge ncbi-datasets-cli
#
# Output:
#   anchors_plastome_universal.fa
#
# References:
#   Li H. 2018. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics.
#   Turudić et al. 2023. Towards the Well-Tempered Chloroplast DNA Sequences.

set -euo pipefail
IFS=$'\n\t'

outdir="${1:-anchors_universal}"
mkdir -p "$outdir"
cd "$outdir"

echo "[INFO] Downloading six universal plastid anchor genes from NC_000932.1 (A. thaliana)..."

# Define associative arrays (requires Bash ≥4)
declare -A start=(
	[psbA]=113761
	[rbcL]=55163
	[matK]=5529
	[ycf1]=131104
	[ndhF]=124796
)

declare -A end=(
	[psbA]=114966
	[rbcL]=56520
	[matK]=7024
	[ycf1]=134431
	[ndhF]=126303
)

# Loop over all gene keys
for g in psbA rbcL matK ycf1 ndhF; do
	echo "[INFO] Fetching $g ..."
	echo efetch -db nucleotide -id NC_000932.1 \
		-seq_start "${start[$g]}" \
		-seq_stop "${end[$g]}" \
		-format fasta ">${g}.fa" 2>&1
	efetch -db nucleotide -id NC_000932.1 \
		-seq_start "${start[$g]}" \
		-seq_stop "${end[$g]}" \
		-format fasta >"${g}.fa"
done

# Combine and clean
# cat psbA.fa trnH-GUG.fa rbcL.fa matK.fa ycf1.fa ndhF.fa >anchors_plastome_universal.fa
cat psbA.fa rbcL.fa matK.fa ycf1.fa ndhF.fa >anchors_plastome_universal.fa

# Validate and summarize
if command -v seqkit >/dev/null 2>&1; then
	seqkit stats anchors_plastome_universal.fa
	seqkit seq -n anchors_plastome_universal.fa >anchors.list
	echo "[OK] Anchors combined → anchors_plastome_universal.fa"
	echo "[OK] Gene list:"
	cat anchors.list
else
	echo "[WARN] seqkit not installed — skipping QC summary."
fi

echo "[DONE] All anchors ready in $PWD/anchors_plastome_universal.fa"
