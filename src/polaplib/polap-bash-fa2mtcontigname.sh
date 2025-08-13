#!/usr/bin/env bash
set -euo pipefail

usage() {
	echo "Usage: $0 [-o output.gfa | -o -] <FASTA>"
	echo "  -o FILE   Write output to FILE"
	echo "  -o -      Write output to stdout"
	exit 1
}

out=""
while getopts ":o:" opt; do
	case "$opt" in
	o) out="$OPTARG" ;;
	*) usage ;;
	esac
done
shift $((OPTIND - 1))

if [[ $# -ne 1 ]]; then
	usage
fi

fa="$1"

_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -n "$out" ]]; then
	awk -f "${_POLAPLIB_DIR}/polap-awk-fa2mtcontigname.awk" "$fa" >"$out"
else
	awk -f "${_POLAPLIB_DIR}/polap-awk-fa2mtcontigname.awk" "$fa"
fi
