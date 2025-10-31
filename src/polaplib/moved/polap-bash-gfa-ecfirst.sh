#!/usr/bin/env bash
# polap-bash-gfa-ecfirst.sh
# Version: v0.2.0
# Ensure ec:i is the FIRST tag on every L-line; optionally drop other tags.
set -euo pipefail
IFS=$'\n\t'

[[ $# -ge 2 ]] || {
	echo "Usage: polap-bash-gfa-ecfirst.sh in.gfa out.gfa [--drop-other-tags]" >&2
	exit 2
}
in="$1"
out="$2"
drop="${3:-}"

command -v python3 >/dev/null 2>&1 || {
	echo "ERR: need python3" >&2
	exit 127
}
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
py="${here}/scripts/polap_py_gfa_ec_first.py"

python3 "$py" --gfa "$in" --out "$out" ${drop:+--drop-other-tags}
echo "[OK] wrote: $out"
