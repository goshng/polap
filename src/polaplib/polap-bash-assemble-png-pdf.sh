#!/usr/bin/env bash
# polap-bash-assemble-png-pdf.sh
# Version : v1.0.0  (2025-10-07)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Build a single-page PDF sheet of genome assembly PNGs from a manifest JSON.
# - PNGs are resolved relative to the manifest’s parent dir (species-root).
# - If a PNG path is blank or missing, the R renderer (make_png_grid.R >= v0.2.6)
#   will draw a red “No assembly” panel (no placeholder dependency).
#
# Requires: bash, Rscript, and (preferably) jq. If jq is absent, we use a safe awk fallback.
#
# Usage:
#   polap-bash-assemble-png-pdf.sh \
#     --manifest md/manifest-v6-some.json \
#     --out md/assemblies-v6-some.pdf \
#     [--organelle mt|pt|both] [--rows N] [--cols N] \
#     [--title "Assemblies (v6)"] \
#     [--page-width-in 8.5] [--page-height-in 11] [--margin-in 0.4] \
#     [--cell-pad 0.03] [--debug]
#
# Notes:
#   - rows/cols are optional hints; the R script guarantees a single page.
#   - --debug keeps temp list files for inspection.

set -euo pipefail
IFS=$'\n\t'

###############################################################################
# Globals & defaults
###############################################################################
_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export _POLAPLIB_DIR

MANIFEST=""
OUTPDF=""
ORG="mt" # mt | pt | both
ROWS=0   # 0 = auto
COLS=0   # 0 = auto
TITLE="Assemblies"
PW=11   # page width (in)
PH=8.5  # page height (in)
MG=0.4  # page margin (in)
CP=0.02 # internal cell padding fraction [0..0.1]
DEBUG=0

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --manifest <manifest.json> --out <out.pdf>
                   [--organelle mt|pt|both] [--rows N] [--cols N]
                   [--title "text"]
                   [--page-width-in IN] [--page-height-in IN] [--margin-in IN]
                   [--cell-pad FRAC] [--debug]
EOF
}

###############################################################################
# Parse arguments
###############################################################################
while (($#)); do
	case "$1" in
	--manifest)
		MANIFEST="${2:?}"
		shift 2
		;;
	--out)
		OUTPDF="${2:?}"
		shift 2
		;;
	--organelle)
		ORG="${2:?}"
		shift 2
		;;
	--rows)
		ROWS="${2:?}"
		shift 2
		;;
	--cols)
		COLS="${2:?}"
		shift 2
		;;
	--title)
		TITLE="${2:?}"
		shift 2
		;;
	--page-width-in)
		PW="${2:?}"
		shift 2
		;;
	--page-height-in)
		PH="${2:?}"
		shift 2
		;;
	--margin-in)
		MG="${2:?}"
		shift 2
		;;
	--cell-pad)
		CP="${2:?}"
		shift 2
		;;
	--debug)
		DEBUG=1
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERR] Unknown option: $1" >&2
		usage
		exit 2
		;;
	esac
done

###############################################################################
# Validate
###############################################################################
[[ -n "$MANIFEST" && -n "$OUTPDF" ]] || {
	echo "[ERR] --manifest and --out are required" >&2
	usage
	exit 2
}
[[ -s "$MANIFEST" ]] || {
	echo "[ERR] no such manifest: $MANIFEST" >&2
	exit 2
}

mkdir -p "$(dirname "$OUTPDF")"

# Species root = parent of manifest dir (…/md -> …/)
SR="$(cd "$(dirname "$MANIFEST")/.." && pwd)"

# Create temp lists
LISTTSV="$(mktemp -t polap-pnglist.XXXXXX.tsv)"
LISTCSV="$(mktemp -t polap-pnglist.XXXXXX.csv)"
if [[ "$DEBUG" -eq 1 ]]; then
	echo "[DEBUG] keeping temp files:"
	echo "  TSV:  $LISTTSV"
	echo "  CSV:  $LISTCSV"
else
	trap 'rm -f "$LISTTSV" "$LISTCSV"' EXIT
fi

echo -e "species\tpng" >"$LISTTSV"

###############################################################################
# Extract species, PNG paths (absolute), one row per item and organelle
###############################################################################
if command -v jq >/dev/null 2>&1; then
	case "$ORG" in
	mt)
		jq -r --arg SR "$SR" '
        .items[] | {species} + {png: (.mt.png // "")}
        | [ .species, (if .png=="" then "" else ($SR + "/" + .png) end) ] | @tsv
      ' "$MANIFEST" >>"$LISTTSV"
		;;
	pt)
		jq -r --arg SR "$SR" '
        .items[] | {species} + {png: (.pt.png // "")}
        | [ .species, (if .png=="" then "" else ($SR + "/" + .png) end) ] | @tsv
      ' "$MANIFEST" >>"$LISTTSV"
		;;
	both)
		jq -r --arg SR "$SR" '
        .items[] as $it |
        [
          [$it.species, ( $it.mt.png // "" )],
          [$it.species, ( $it.pt.png // "" )]
        ] | .[]
        | [ .[0], (if .[1]=="" then "" else ($SR + "/" + .[1]) end) ] | @tsv
      ' "$MANIFEST" >>"$LISTTSV"
		;;
	*)
		echo "[ERR] --organelle must be mt|pt|both" >&2
		exit 2
		;;
	esac
else
	# awk fallback (split-by-quote, no fragile regex escapes)
	awk -F'"' -v ORG="$ORG" -v SR="$SR" '
    /"species"[[:space:]]*:/ { species=$4; mt=""; pt=""; start=1 }
    /"png"[[:space:]]*:/ {
      p=$4; if (p !~ /^\//) p=SR "/" p
      if (p ~ /\/mt[._]/) mt=p
      if (p ~ /\/pt[._]/) pt=p
    }
    /},[[:space:]]*{/ || /\}[[:space:]]*\][[:space:]]*\}/ {
      if (start && species!="") {
        if (ORG=="mt")      printf "%s\t%s\n", species, mt
        else if (ORG=="pt") printf "%s\t%s\n", species, pt
        else { printf "%s\t%s\n%s\t%s\n", species, mt, species, pt }
      }
      start=0
    }
  ' "$MANIFEST" >>"$LISTTSV"
fi

# TSV -> CSV (keep empty png fields; R handles “No assembly”)
awk -F'\t' 'BEGIN{OFS=","} NR==1{print "species","png"; next} {print $1,$2}' "$LISTTSV" >"$LISTCSV"

echo "[INFO] First few PNG entries:"
sed -n '1,8p' "$LISTCSV"

###############################################################################
# Render (R script >= v0.2.6 handles margins, single page, and “No assembly”)
###############################################################################
Rscript "${_POLAPLIB_DIR}/scripts/make_png_grid.R" \
	--list "$LISTCSV" \
	--out "$OUTPDF" \
	--title "$TITLE" \
	--page-width-in "$PW" \
	--page-height-in "$PH" \
	--margin-in "$MG" \
	--gap-mm 5 \
	--rows "$ROWS" \
	--label-cex 0.5 \
	--cols "$COLS"

echo "[OK] Wrote: $OUTPDF"
