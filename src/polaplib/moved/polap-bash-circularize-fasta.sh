#!/usr/bin/env bash
# polap-bash-circularize-fasta.sh
# Version: v0.2.0
# Turn a single-record FASTA (linear path) into a circular sequence by trimming head/tail exact overlap.
set -euo pipefail
IFS=$'\n\t'

in="" out="" min_ovl=2000 max_ovl=20000 rotate=""

usage() {
	cat <<'USAGE'
Usage:
  polap-bash-circularize-fasta.sh --in in.fa --out out.circular.fa [--min-ovl 2000] [--max-ovl 20000] [--rotate-seed SEQ]
USAGE
}

while (($#)); do
	case "$1" in
	--in)
		in="$2"
		shift 2
		;;
	--out)
		out="$2"
		shift 2
		;;
	--min-ovl)
		min_ovl="$2"
		shift 2
		;;
	--max-ovl)
		max_ovl="$2"
		shift 2
		;;
	--rotate-seed)
		rotate="$2"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "ERR: unknown arg $1" >&2
		exit 2
		;;
	esac
done

[[ -s "$in" && -n "$out" ]] || {
	echo "ERR: need --in and --out" >&2
	exit 2
}
command -v python3 >/dev/null || {
	echo "ERR: need python3" >&2
	exit 127
}

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
py="${here}/scripts/circularize_fasta_exact.py"

python3 "$py" --in "$in" --out "$out" --min-ovl "$min_ovl" --max-ovl "$max_ovl" ${rotate:+--rotate-seed "$rotate"}
echo "[OK] wrote: $out"
