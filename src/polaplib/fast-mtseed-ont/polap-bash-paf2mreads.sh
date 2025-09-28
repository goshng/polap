#!/usr/bin/env bash
# polap-bash-paf2mreads.sh  v0.0.1
# Usage: paf awk_script seedreads_py out_ids
set -euo pipefail
paf="$1"; awk_script="$2"; seed_py="$3"; out="$4"
awk -f "$awk_script" "$paf" | python "$seed_py" > "$out"
