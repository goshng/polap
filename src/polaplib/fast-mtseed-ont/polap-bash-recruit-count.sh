#!/usr/bin/env bash
# polap-bash-recruit-count.sh  v0.0.1
# Usage: polap-bash-recruit-count.sh in.bam out.tsv
set -euo pipefail
bam="$1"; out="$2"
samtools view -F 2304 "$bam" \
  | awk '{print $3, $1}' \
  | sort -u \
  | awk '{c[$1]++} END{for(k in c) print k"\t"c[k]}' \
  | sort -k2,2nr > "$out"
