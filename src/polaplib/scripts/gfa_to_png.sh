#!/usr/bin/env bash
################################################################################
# scripts/gfa_to_png-v0.1.0.sh
#
# Version : v0.1.0
# Author  : Sang Chul Choi (POLAP)
# Date    : 2025-10-05
# License : GPL-3.0+
#
# Purpose :
#   Render assembly GFA to PNG via Bandage through polap.sh, with caching and a
#   graceful fallback to placeholder when Bandage is unavailable.
#
# Usage   :
#   gfa_to_png-v0.1.0.sh <species> [pt|mt] [bandage|no-bandage]
################################################################################
set -euo pipefail
IFS=$'\n\t'

_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
export _POLAPLIB_DIR

sp="${1:?species}"
kind="${2:-pt}"     # pt|mt
how="${3:-bandage}" # bandage|no-bandage

t="t5"
inum="0"
base="${sp}/${t}/${inum}"

if [[ "$kind" == "pt" ]]; then
	gfa="${base}/polap-readassemble-1-pt/pt.1.gfa"
	png="${base}/polap-readassemble-1-pt/pt.1.png"
else
	gfa="${base}/polap-readassemble-1-mt/mt.1.gfa"
	png="${base}/polap-readassemble-1-mt/mt.1.png"
fi

mkdir -p "$(dirname "$png")"
# Cache: skip if png newer than gfa
if [[ -s "$png" && -s "$gfa" && "$png" -nt "$gfa" ]]; then
	exit 0
fi

if [[ -s "$gfa" ]]; then
	if [[ "$how" == "bandage" ]] && command -v polap.sh >/dev/null 2>&1; then
		if ! polap.sh bandage png "$gfa" "$png"; then
			cp -f "${_POLAPLIB_DIR}/man/figures/na.png" "$png" 2>/dev/null || true
		fi
	else
		# placeholder (requires ImageMagick convert if present)
		if command -v convert >/dev/null 2>&1; then
			convert -size 800x400 xc:white -gravity center -pointsize 18 \
				-annotate 0 "${sp} (${kind})" "$png" 2>/dev/null ||
				cp -f "${_POLAPLIB_DIR}/man/figures/na.png" "$png" 2>/dev/null || true
		else
			cp -f "${_POLAPLIB_DIR}/man/figures/na.png" "$png" 2>/dev/null || true
		fi
	fi
else
	cp -f "${_POLAPLIB_DIR}/man/figures/na.png" "$png" 2>/dev/null || true
fi
