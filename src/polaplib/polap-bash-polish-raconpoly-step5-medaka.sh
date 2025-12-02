#!/usr/bin/env bash
# File: polap-bash-polish-raconpoly-step5-medaka.sh
# Version: v0.9.2
set -Eeuo pipefail
# Args: -R reads.fq.gz -A asm.fa -m model -t threads -o outdir
while getopts ":R:A:m:t:o:" opt; do
	case "$opt" in
	R) READS="$OPTARG" ;;
	A) ASM="$OPTARG" ;;
	m) MODEL="$OPTARG" ;;
	t) THREADS="$OPTARG" ;;
	o) OUTDIR="$OPTARG" ;;
	*)
		echo "bad opt"
		exit 2
		;;
	esac
done
[[ -s "$READS" && -s "$ASM" && -n "$OUTDIR" ]] || {
	echo "args?"
	exit 2
}
install -d -- "$OUTDIR"
minimap2 -ax map-ont -t "$THREADS" "$ASM" "$READS" | samtools sort -@ "$THREADS" -o "$OUTDIR/reads.sorted.bam" -
samtools index -@ "$THREADS" "$OUTDIR/reads.sorted.bam"
if ! medaka_consensus -i "$READS" -d "$ASM" -o "$OUTDIR" -m "$MODEL" -t "$THREADS"; then
	echo "medaka_consensus failed; continuing without medaka" >&2
fi
