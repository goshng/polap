#!/usr/bin/env bash
# File: polap-bash-polish-raconpoly-step8-polypolish.sh
# Version: v0.9.2
set -Eeuo pipefail
# Args: -a aligner -S mm2_secondaries(0|1) -A asm.fa -1 R1 -2 R2 -B allow.bed -t threads -d pp_dir -o out.fa -b out.name.bam
while getopts ":a:S:A:1:2:B:t:d:o:b:" opt; do
	case "$opt" in
	a) ALIGNER="$OPTARG" ;;
	S) MM2SEC="$OPTARG" ;;
	A) ASM="$OPTARG" ;;
	1) R1="$OPTARG" ;;
	2) R2="$OPTARG" ;;
	B) ALLOW="$OPTARG" ;;
	t) THREADS="$OPTARG" ;;
	d) PPD="$OPTARG" ;;
	o) OUTFA="$OPTARG" ;;
	b) OUTBAM="$OPTARG" ;;
	*)
		echo "bad opt"
		exit 2
		;;
	esac
done
[[ -s "$ASM" && -s "$R1" && -s "$R2" && -n "$OUTFA" && -n "$OUTBAM" && -n "$PPD" ]] || {
	echo "args?"
	exit 2
}
install -d -- "$PPD/bt2idx"

case "$ALIGNER" in
bowtie2)
	bowtie2-build "$ASM" "$PPD/bt2idx/asm" >"$PPD/pp_build.out" 2>"$PPD/pp_build.err"
	bowtie2 --very-sensitive -a -x "$PPD/bt2idx/asm" -1 "$R1" -2 "$R2" -p "$THREADS" 2>"$PPD/pp_bt2.err" |
		samtools view -bh - |
		samtools view -bh -L "$ALLOW" - |
		samtools sort -n -@ "$THREADS" -o "$OUTBAM" - 2>"$PPD/pp_sortn.err"
	;;
bwa)
	bwa-mem2 index "$ASM" >"$PPD/pp_bwa_index.out" 2>"$PPD/pp_bwa_index.err"
	bwa-mem2 mem -t "$THREADS" -a "$ASM" "$R1" "$R2" 2>"$PPD/pp_bwa.err" |
		samtools view -bh - |
		samtools view -bh -L "$ALLOW" - |
		samtools sort -n -@ "$THREADS" -o "$OUTBAM" - 2>"$PPD/pp_sortn.err"
	;;
minimap2)
	mm2_extra=""
	[[ "$MM2SEC" == "1" ]] && mm2_extra="--secondary=yes -N 50"
	minimap2 -x sr -a $mm2_extra -t "$THREADS" "$ASM" "$R1" "$R2" 2>"$PPD/pp_mmsr.err" |
		samtools view -bh - |
		samtools view -bh -L "$ALLOW" - |
		samtools sort -n -@ "$THREADS" -o "$OUTBAM" - 2>"$PPD/pp_sortn.err"
	;;
*)
	echo "bad aligner"
	exit 2
	;;
esac

samtools view "$OUTBAM" | polypolish polish "$ASM" /dev/stdin >"$OUTFA" 2>"$PPD/polypolish.err" || true
