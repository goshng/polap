#!/usr/bin/env bash
# polap-bash-gfa-fix-ec-first.sh
# Version: v0.1.0
# Make sure every L line has ec:i AND that ec:i is the FIRST tag after the CIGAR.
# Optionally drop all other edge tags (safer for gfatk 0.2.2).

set -euo pipefail
IFS=$'\n\t'

if (($# < 2)); then
	echo "Usage: polap-bash-gfa-fix-ec-first.sh in.gfa out.gfa [--drop-other-tags]" >&2
	exit 2
fi

in="$1"
out="$2"
drop="${3:-}"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
py="${script_dir}/scripts/polap_py_gfa_ec_first.py"

python3 "$py" --gfa "$in" --out "$out" ${drop:+--drop-other-tags}
echo "[OK] Wrote: $out"
