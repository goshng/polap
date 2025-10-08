#!/usr/bin/env bash
################################################################################
# scripts/build_manifest-v0.1.0.sh
#
# Version : v0.1.0
# Author  : Sang Chul Choi (POLAP)
# Date    : 2025-10-05
# License : GPL-3.0+
#
# Purpose :
#   Emit a simple manifest JSON of generated figure PNGs for a given set/type.
#
# Usage   :
#   build_manifest-v0.1.0.sh --set SET --type pt|mt --out DIR
################################################################################
set -euo pipefail
IFS=$'\n\t'

_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
export _POLAPLIB_DIR

SET="some"
TYPE="pt"
OUT="${_POLAPLIB_DIR}/man/v0.5.4/figures"

usage() {
	sed -n '1,80p' "${BASH_SOURCE[0]}"
}

while (($#)); do
	case "$1" in
	--set)
		SET="${2:?}"
		shift 2
		;;
	--type)
		TYPE="${2:?}"
		shift 2
		;;
	--out)
		OUT="${2:?}"
		shift 2
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERR] Unknown option: $1" >&2
		exit 2
		;;
	esac
done

mkdir -p "${OUT}"
manifest="${OUT}/manifest-${SET}-${TYPE}.json"
tmp="$(mktemp -t polap-manifest.XXXXXX)"

{
	printf '{\n  "set":"%s", "type":"%s", "figures":[\n' "${SET}" "${TYPE}"
	first=1
	# We can only rely on directory structure; read species from CSV if present.
	if [[ -s "${_POLAPLIB_DIR}/polap-bash-report.csv" ]]; then
		species_list=$(csvtk -t cut -f species "${_POLAPLIB_DIR}/polap-bash-report.csv" | tail -n +2 2>/dev/null || true)
	else
		species_list=""
	fi
	if [[ -n "${species_list}" ]]; then
		while read -r sp; do
			base="${sp}/t5/0"
			if [[ "${TYPE}" == "pt" ]]; then
				png="${base}/polap-readassemble-1-pt/pt.1.png"
			else
				png="${base}/polap-readassemble-1-mt/mt.1.png"
			fi
			[[ $first -eq 1 ]] || printf ",\n"
			first=0
			printf '    {"species":"%s","png":"%s"}' "$sp" "$png"
		done <<<"${species_list}"
	fi
	printf '\n  ]\n}\n'
} >"${tmp}"

mv "${tmp}" "${manifest}"
echo "${manifest}"
