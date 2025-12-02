#!/usr/bin/env bash
# File: polap-bash-polish-raconpoly-step4-racon-round.sh
# Version: v0.9.2
set -Eeuo pipefail
# Args: -R reads.fq.gz -A asm.fa -p PRESET -t THREADS -i MINID -a MINALEN -n ROUND -O out.fa -D stage_dir
while getopts ":R:A:p:t:i:a:n:O:D:" opt; do
	case "$opt" in
	R) READS="$OPTARG" ;;
	A) ASM="$OPTARG" ;;
	p) PRESET="$OPTARG" ;;
	t) THREADS="$OPTARG" ;;
	i) MINID="$OPTARG" ;;
	a) MINALEN="$OPTARG" ;;
	n) ROUND="$OPTARG" ;;
	O) OUTFA="$OPTARG" ;;
	D) STAGED="$OPTARG" ;;
	*)
		echo "bad opt"
		exit 2
		;;
	esac
done
[[ -s "$READS" && -s "$ASM" && -n "$OUTFA" && -n "$STAGED" && -n "$ROUND" ]] || {
	echo "args?"
	exit 2
}
POLAPLIB_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
SCRIPTS_DIR="${POLAPLIB_DIR}/scripts"
PAF_FILTER="${SCRIPTS_DIR}/polap_py_paf_filter.py"
RAW="$STAGED/r${ROUND}.raw.paf"
FLT="$STAGED/r${ROUND}.flt.paf"
minimap2 -x "$PRESET" -t "$THREADS" -c --secondary=no "$ASM" "$READS" >"$RAW" 2>"$STAGED/r${ROUND}.minimap2.err"
python3 "$PAF_FILTER" --min-ident "$MINID" --min-alen "$MINALEN" --in "$RAW" --out "$FLT" 2>>"$STAGED/r${ROUND}.paf_filter.err"
[[ -s "$FLT" ]] || cp -f "$RAW" "$FLT"
racon -t "$THREADS" "$READS" "$FLT" "$ASM" >"$OUTFA" 2>"$STAGED/r${ROUND}.racon.err"
