#!/usr/bin/env bash
# File: polap-bash-polish-raconpoly-step1-assign-ont.sh
# Version: v0.9.2
set -Eeuo pipefail
# Args: -r ref.fa -q ont.fq.gz -x target_prefix -m MAPQ -p ONT_PRESET -t THREADS -o out_ids
while getopts ":r:q:x:m:p:t:o:" opt; do
	case "$opt" in
	r) REF="$OPTARG" ;;
	q) FQ="$OPTARG" ;;
	x) TGT_PREFIX="$OPTARG" ;;
	m) MAPQ="$OPTARG" ;;
	p) PRESET="$OPTARG" ;;
	t) THREADS="$OPTARG" ;;
	o) OUT="$OPTARG" ;;
	*)
		echo "bad opt"
		exit 2
		;;
	esac
done
[[ -s "$REF" && -s "$FQ" && -n "$TGT_PREFIX" && -n "$OUT" ]] || {
	echo "args?"
	exit 2
}

POLAPLIB_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
SCRIPTS_DIR="${POLAPLIB_DIR}/scripts"
ASSIGN_ONT_IDS="${SCRIPTS_DIR}/polap_py_assign_ont_ids.py"
minimap2 -ax "$PRESET" -t "$THREADS" "$REF" "$FQ" |
	samtools view -h -F 0x900 - |
	python3 "$ASSIGN_ONT_IDS" --sam - --target-prefix "$TGT_PREFIX" --mapq "$MAPQ" \
		>"$OUT"
