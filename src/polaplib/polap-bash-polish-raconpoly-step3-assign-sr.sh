#!/usr/bin/env bash
# File: polap-bash-polish-raconpoly-step3-assign-sr.sh
# Version: v0.9.2
set -Eeuo pipefail
# Args: -r combined.fa -x target_prefix -m MAPQ -t THREADS -1 R1 -2 R2 -o out_ids
while getopts ":r:x:m:t:1:2:o:" opt; do
	case "$opt" in
	r) REF="$OPTARG" ;;
	x) TGT_PREFIX="$OPTARG" ;;
	m) MAPQ="$OPTARG" ;;
	t) THREADS="$OPTARG" ;;
	1) R1="$OPTARG" ;;
	2) R2="$OPTARG" ;;
	o) OUT="$OPTARG" ;;
	*)
		echo "bad opt"
		exit 2
		;;
	esac
done
[[ -s "$REF" && -s "$R1" && -s "$R2" && -n "$OUT" ]] || {
	echo "args?"
	exit 2
}
POLAPLIB_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
SCRIPTS_DIR="${POLAPLIB_DIR}/scripts"
ASSIGN_SR_IDS="${SCRIPTS_DIR}/polap_py_assign_sr_ids.py"
IDXDIR="$(dirname -- "$OUT")/bt2idx"
install -d -- "$IDXDIR"
bowtie2-build "$REF" "$IDXDIR/comb"
bowtie2 --very-sensitive -x "$IDXDIR/comb" -1 "$R1" -2 "$R2" -p "$THREADS" |
	samtools view -bh - |
	python3 "$ASSIGN_SR_IDS" --bam - --target-prefix "$TGT_PREFIX" --mapq "$MAPQ" \
		>"$OUT"
