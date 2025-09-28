#!/usr/bin/env bash
set -euo pipefail

# polap-bash-reads.sh — fast organelle-read table generator (C first, R fallback)
# Usage:
#   polap-bash-reads.sh --mt mt.paf --pt pt.paf --output outdir \
#     --min-mapq 1 --min-pt 0 --min-identity 0.15

usage() {
	cat <<EOF
Usage: $0 --mt FILE --pt FILE --output DIR [--min-mapq INT] [--min-pt INT] [--min-identity FLOAT]
EOF
}

# ---------- parse args ----------
MT=""
PT=""
OUT=""
MIN_MAPQ=1
MIN_PT=0
MIN_ID=0.15

while [[ $# -gt 0 ]]; do
	case "$1" in
	--mt)
		MT="$2"
		shift 2
		;;
	--pt)
		PT="$2"
		shift 2
		;;
	--output)
		OUT="$2"
		shift 2
		;;
	--min-mapq)
		MIN_MAPQ="$2"
		shift 2
		;;
	--min-pt)
		MIN_PT="$2"
		shift 2
		;;
	--min-identity)
		MIN_ID="$2"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[err] Unknown arg: $1"
		usage
		exit 1
		;;
	esac
done

[[ -z "$MT" || -z "$PT" || -z "$OUT" ]] && {
	usage
	exit 1
}
mkdir -p "$OUT"

# ---------- try C version ----------
bin="./polap-reads"
src_dir="$(dirname "${BASH_SOURCE[0]}")/scripts"
c_src="${src_dir}/polap-c-reads.c"

if [[ ! -x "$bin" ]]; then
	if command -v cc >/dev/null 2>&1; then
		echo "[info] building C binary: $bin"
		cc -O3 -march=native -pipe -o "$bin" "$c_src"
	fi
fi

if [[ -x "$bin" ]]; then
	"$bin" \
		--mt "$MT" --pt "$PT" --output "$OUT" \
		--min-mapq "$MIN_MAPQ" --min-pt "$MIN_PT" --min-identity "$MIN_ID"
	exit $?
fi

# ---------- fallback to fast R ----------
if command -v Rscript >/dev/null 2>&1; then
	Rscript "${src_dir}/polap-r-reads-fast.R" \
		--mt "$MT" --pt "$PT" --output "$OUT" \
		--min-mapq "$MIN_MAPQ" --min-pt "$MIN_PT" --min-identity "$MIN_ID"
	exit $?
else
	echo "[err] Neither compiled C binary nor Rscript available."
	exit 1
fi
