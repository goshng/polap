#!/usr/bin/env bash
set -euo pipefail

# estimate_genome_size_cov.sh
#
# Coverage-based genome size estimate from:
#   (A) mosdepth summary (depth.summary.txt)
#   (B) the FASTQ reads used for mapping (e.g., >=1kb filtered)
#
# Requires helper scripts placed next to this file:
#   compute_weighted_median.awk
#   flag_organelle.awk
#   sum_len_times_mean.awk
#   filter_nuclear_only.awk
#   fastq_total_bases.awk
#   Gcov.py
#
# Usage:
#   ./estimate_genome_size_cov.sh \
#     -d depth.summary.txt \
#     -r clr.min1k.fq \
#     -o outdir \
#     [-m 10000] \
#     [-M 2.5]
#
# Options:
#   -d  Path to mosdepth summary (TSV with header: chrom length bases mean min max)
#   -r  Path to FASTQ used for mapping (the same reads used to make the BAM)
#   -o  Output directory (all outputs go here)
#   -m  MINLEN: min contig length to consider for medians (default: 10000)
#   -M  ORG_MULT: organelle/high-copy depth multiplier (default: 2.5)
#
# Outputs (all under -o):
#   pivot1.txt
#   summary.with_flag.tsv
#   nuclear.summary.tsv
#   pivot_nuclear.txt
#   mapped_all.txt
#   mapped_nuc.txt
#   total_bases.txt
#   Gcov.txt              # final coverage-based genome size (bp)
#   estimate.summary.tsv  # small human-readable summary
#
# Prints a short summary to stdout as well.

DEPTH_SUMMARY=""
READS_FASTQ=""
OUTDIR=""
MINLEN=10000
ORG_MULT=2.5

while getopts ":d:r:o:m:M:h" opt; do
	case "$opt" in
	d) DEPTH_SUMMARY="$OPTARG" ;;
	r) READS_FASTQ="$OPTARG" ;;
	o) OUTDIR="$OPTARG" ;;
	m) MINLEN="$OPTARG" ;;
	M) ORG_MULT="$OPTARG" ;;
	h)
		sed -n '1,200p' "$0" | sed -n '1,80p'
		exit 0
		;;
	\?)
		echo "ERROR: Invalid option -$OPTARG" >&2
		exit 2
		;;
	:)
		echo "ERROR: Option -$OPTARG requires an argument" >&2
		exit 2
		;;
	esac
done

# Basic checks
if [[ -z "$DEPTH_SUMMARY" || -z "$READS_FASTQ" || -z "$OUTDIR" ]]; then
	echo "ERROR: -d depth.summary.txt, -r reads.fastq, and -o outdir are required." >&2
	exit 2
fi
if [[ ! -s "$DEPTH_SUMMARY" ]]; then
	echo "ERROR: Depth summary not found or empty: $DEPTH_SUMMARY" >&2
	exit 2
fi
if [[ ! -s "$READS_FASTQ" ]]; then
	echo "ERROR: FASTQ not found or empty: $READS_FASTQ" >&2
	exit 2
fi

mkdir -p "$OUTDIR"

# Resolve helper scripts relative to this script's directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
AWK_MED="$SCRIPT_DIR/compute_weighted_median.awk"
AWK_FLAG="$SCRIPT_DIR/flag_organelle.awk"
AWK_SUM="$SCRIPT_DIR/sum_len_times_mean.awk"
AWK_NUC="$SCRIPT_DIR/filter_nuclear_only.awk"
AWK_FQ="$SCRIPT_DIR/fastq_total_bases.awk"
PY_GCOV="$SCRIPT_DIR/Gcov.py"

for f in "$AWK_MED" "$AWK_FLAG" "$AWK_SUM" "$AWK_NUC" "$AWK_FQ" "$PY_GCOV"; do
	[[ -s "$f" ]] || {
		echo "ERROR: Helper missing: $f" >&2
		exit 2
	}
done

# A) preliminary pivot (weighted median depth across contigs ≥ MINLEN)
sort -k4,4n "$DEPTH_SUMMARY" |
	awk -v minlen="$MINLEN" -f "$AWK_MED" \
		>"$OUTDIR/pivot1.txt"
pivot1="$(cat "$OUTDIR/pivot1.txt")"

# B) flag organelle/high-copy & write summary.with_flag.tsv
awk -v pivot="$pivot1" -v mult="$ORG_MULT" -v minlen="$MINLEN" \
	-f "$AWK_FLAG" "$DEPTH_SUMMARY" \
	>"$OUTDIR/summary.with_flag.tsv"

# C) nuclear-only subset
awk -v minlen="$MINLEN" -f "$AWK_NUC" "$OUTDIR/summary.with_flag.tsv" \
	>"$OUTDIR/nuclear.summary.tsv"

# D) recompute nuclear pivot (weighted median) on nuclear subset
sort -k4,4n "$OUTDIR/nuclear.summary.tsv" |
	awk -v minlen="$MINLEN" -f "$AWK_MED" \
		>"$OUTDIR/pivot_nuclear.txt"
Dn="$(cat "$OUTDIR/pivot_nuclear.txt")"

# E) mapped bases: sum(length*mean) for all vs nuclear-only
awk -f "$AWK_SUM" "$DEPTH_SUMMARY" >"$OUTDIR/mapped_all.txt"
awk -f "$AWK_SUM" "$OUTDIR/nuclear.summary.tsv" >"$OUTDIR/mapped_nuc.txt"
frac_nuc="$(
	awk 'NR==FNR{a=$1; next} {printf("%.6f\n", (a>0 ? $1/a : 0))}' \
		"$OUTDIR/mapped_all.txt" "$OUTDIR/mapped_nuc.txt"
)"

# F) total read bases used (FASTQ)
awk -f "$AWK_FQ" "$READS_FASTQ" >"$OUTDIR/total_bases.txt"
Tb="$(cat "$OUTDIR/total_bases.txt")"

# G) genome size estimate
python "$PY_GCOV" "$OUTDIR/total_bases.txt" "$OUTDIR/pivot_nuclear.txt" "$frac_nuc" \
	>"$OUTDIR/Gcov.txt"

# Summary TSV
{
	echo -e "metric\tvalue"
	echo -e "MINLEN\t$MINLEN"
	echo -e "ORG_MULT\t$ORG_MULT"
	echo -e "pivot1_depth\t$pivot1"
	echo -e "nuclear_depth_Dn\t$Dn"
	echo -e "frac_nuclear\t$frac_nuc"
	echo -e "total_bases\t$Tb"
	echo -e "Gcov_bp\t$(cat "$OUTDIR/Gcov.txt")"
} >"$OUTDIR/estimate.summary.tsv"

# Console summary
echo "=== Coverage-based genome size estimate ==="
echo "  Output dir        : $OUTDIR"
echo "  MINLEN            : $MINLEN"
echo "  ORG_MULT          : $ORG_MULT"
echo "  pivot1 depth      : $pivot1"
echo "  nuclear depth Dn  : $Dn"
echo "  frac_nuclear      : $frac_nuc"
echo "  total bases (Tb)  : $Tb"
echo "  Gcov (bp)         : $(cat "$OUTDIR/Gcov.txt")"
echo "  Summary TSV       : $OUTDIR/estimate.summary.tsv"
