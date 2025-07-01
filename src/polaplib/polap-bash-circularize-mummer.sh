#!/bin/bash

set -e

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.fasta> [min_length=1000] [min_identity=95]"
  exit 1
fi

FA="$1"
MINLEN="${2:-1000}" # Minimum overlap length (default 1000)
MINID="${3:-95}"    # Minimum percent identity (default 95)
PREFIX="circularize_tmp"
OUT="${FA%.fasta}.circularized.fasta"

echo "[INFO] Input: $FA"
echo "[INFO] Searching for imperfect terminal direct repeats..."

# 1. Self-alignment using MUMmer
nucmer --maxmatch -p "$PREFIX" "$FA" "$FA" >/dev/null
show-coords -rcl "$PREFIX.delta" >"$PREFIX.coords"

# 2. Get sequence length
SEQLEN=$(seqkit stats "$FA" | awk 'NR==2 {print $5}')

# 3. Filter for terminal-to-terminal matches
REPEAT=$(awk -v minlen="$MINLEN" -v minid="$MINID" -v seqlen="$SEQLEN" 'NR > 5 && $7 >= minid && $5 >= minlen && (($1 <= minlen && $2 >= seqlen - minlen) || ($2 <= minlen && $1 >= seqlen - minlen)) {print $1, $2, $5, $7; exit}' "$PREFIX.coords")

if [[ -z "$REPEAT" ]]; then
  echo "[WARN] No suitable head–tail repeat found (≥$MINLEN bp, ≥$MINID% ID)"
  exit 1
fi

# 4. Parse and trim
read START END LEN PID <<<"$REPEAT"
TRIMLEN="$LEN"
NEWLEN=$((SEQLEN - TRIMLEN))

echo "[INFO] Found terminal match: $TRIMLEN bp at $START–$END, $PID% identity"
echo "[INFO] Circularizing by trimming one copy (keeping 1..$NEWLEN bp)"

# 5. Output circularized sequence
seqkit subseq --region 1:"$NEWLEN" "$FA" >"$OUT"

echo "[DONE] Saved circularized sequence to $OUT"

# Cleanup
# rm -f "$PREFIX."*
