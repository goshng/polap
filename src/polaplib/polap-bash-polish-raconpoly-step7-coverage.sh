#!/usr/bin/env bash
# File: polap-bash-polish-raconpoly-step7-coverage.sh
# Version: v0.9.2
set -Eeuo pipefail
# Args: -A asm.fa -Q allONT.fq -E assignedONT.fq -t threads -b bin -q minMAPQ -m mask.bed -o cov_dir -R cov_plot.R
while getopts ":A:Q:E:t:b:q:m:o:R:" opt; do
	case "$opt" in
	A) ASM="$OPTARG" ;;
	Q) ALL="$OPTARG" ;;
	E) ASN="$OPTARG" ;;
	t) THREADS="$OPTARG" ;;
	b) BIN="$OPTARG" ;;
	q) MINQ="$OPTARG" ;;
	m) MASK="$OPTARG" ;;
	o) COVD="$OPTARG" ;;
	R) COVPLOT="$OPTARG" ;;
	*)
		echo "bad opt"
		exit 2
		;;
	esac
done
[[ -s "$ASM" && -s "$ALL" && -s "$ASN" && -n "$COVD" ]] || {
	echo "args?"
	exit 2
}
install -d -- "$COVD"
POLAPLIB_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
SCRIPTS_DIR="${POLAPLIB_DIR}/scripts"
COV_BINS="${SCRIPTS_DIR}/polap_py_bam_coverage_bins.py"

minimap2 -ax map-ont -t "$THREADS" "$ASM" "$ALL" | samtools sort -@ "$THREADS" -o "$COVD/allONT.sorted.bam" -
samtools index -@ "$THREADS" "$COVD/allONT.sorted.bam"
minimap2 -ax map-ont -t "$THREADS" "$ASM" "$ASN" | samtools sort -@ "$THREADS" -o "$COVD/assignedONT.sorted.bam" -
samtools index -@ "$THREADS" "$COVD/assignedONT.sorted.bam"

python3 "$COV_BINS" --fasta "$ASM" --bam "$COVD/allONT.sorted.bam" \
	--bin "$BIN" --min-mapq "$MINQ" --source allONT \
	--out-tsv "$COVD/coverage.allONT.bin${BIN}.tsv"
python3 "$COV_BINS" --fasta "$ASM" --bam "$COVD/assignedONT.sorted.bam" \
	--bin "$BIN" --min-mapq "$MINQ" --source assignedONT \
	--out-tsv "$COVD/coverage.assignedONT.bin${BIN}.tsv"

if command -v Rscript >/dev/null 2>&1 && [[ -s "$COVPLOT" ]]; then
	Rscript --vanilla "$COVPLOT" \
		--fasta "$ASM" \
		--cov-all "$COVD/coverage.allONT.bin${BIN}.tsv" \
		--cov-assigned "$COVD/coverage.assignedONT.bin${BIN}.tsv" \
		--mask "$MASK" \
		--bin "$BIN" \
		--out-pdf "$COVD/coverage_overlay.pdf" \
		--out-tsv "$COVD/coverage_summary.tsv" \
		--title "coverage (all ONT vs assigned ONT)"
fi
