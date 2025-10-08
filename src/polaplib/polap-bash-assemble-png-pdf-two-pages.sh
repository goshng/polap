#!/usr/bin/env bash
# polap-bash-assemble-png-pdf-two-pages.sh
# Version : v1.0.0  (2025-10-07)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Build a two-page PDF: page 1 = pt, page 2 = mt, from a single manifest JSON.
# PNGs are resolved relative to the manifest’s parent dir (species-root).
# Missing PNGs are left empty; the R renderer draws red "No assembly" panels.

set -euo pipefail
IFS=$'\n\t'

_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export _POLAPLIB_DIR

# ---------------- defaults ----------------
MANIFEST=""
OUTPDF=""
ROWS=0
COLS=0
PW=11
PH=8.5
MG=0.45
GAP=4.0
LABEL_CEX=0.60
TITLE_CEX=1.0
LABEL_STRIP=0.16
PT_TITLE="Plastid assemblies (pt)"
MT_TITLE="Mitochondrial assemblies (mt)"
DEBUG=0

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --manifest <manifest.json> --out <out.pdf>
                   [--rows N] [--cols N]
                   [--page-width-in IN] [--page-height-in IN] [--margin-in IN]
                   [--gap-mm MM] [--label-cex CEX] [--title-cex CEX]
                   [--label-strip-frac FRAC]
                   [--pt-title "text"] [--mt-title "text"]
                   [--debug]
EOF
}

# --------------- parse args ---------------
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
	--rows)
		ROWS="${2:?}"
		shift 2
		;;
	--cols)
		COLS="${2:?}"
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
	--gap-mm)
		GAP="${2:?}"
		shift 2
		;;
	--label-cex)
		LABEL_CEX="${2:?}"
		shift 2
		;;
	--title-cex)
		TITLE_CEX="${2:?}"
		shift 2
		;;
	--label-strip-frac)
		LABEL_STRIP="${2:?}"
		shift 2
		;;
	--pt-title)
		PT_TITLE="${2:?}"
		shift 2
		;;
	--mt-title)
		MT_TITLE="${2:?}"
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

# Species root = parent of manifest dir (…/md → …/)
SR="$(cd "$(dirname "$MANIFEST")/.." && pwd)"

# temp lists
PT_TSV="$(mktemp -t polap-pt.XXXXXX.tsv)"
MT_TSV="$(mktemp -t polap-mt.XXXXXX.tsv)"
PT_CSV="$(mktemp -t polap-pt.XXXXXX.csv)"
MT_CSV="$(mktemp -t polap-mt.XXXXXX.csv)"
if [[ "$DEBUG" -eq 1 ]]; then
	echo "[DEBUG] keeping temp files:"
	echo "  $PT_TSV"
	echo "  $MT_TSV"
	echo "  $PT_CSV"
	echo "  $MT_CSV"
else
	trap 'rm -f "$PT_TSV" "$MT_TSV" "$PT_CSV" "$MT_CSV"' EXIT
fi

echo -e "species\tpng" >"$PT_TSV"
echo -e "species\tpng" >"$MT_TSV"

# --------- extract lists (include empties so all species appear) ----------
if command -v jq >/dev/null 2>&1; then
	jq -r --arg SR "$SR" '
    .items[] | {species} + {png: (.pt.png // "")}
    | [ .species, (if .png=="" then "" else ($SR + "/" + .png) end) ] | @tsv
  ' "$MANIFEST" >>"$PT_TSV"

	jq -r --arg SR "$SR" '
    .items[] | {species} + {png: (.mt.png // "")}
    | [ .species, (if .png=="" then "" else ($SR + "/" + .png) end) ] | @tsv
  ' "$MANIFEST" >>"$MT_TSV"

else
	# awk fallback (split by quotes; record each block’s pt and mt)
	awk -F'"' -v SR="$SR" '
    /"species"[[:space:]]*:/ { species=$4; mt=""; pt=""; start=1 }
    /"png"[[:space:]]*:/ {
      p=$4; if (p !~ /^\//) p=SR "/" p
      if (p ~ /\/mt[._]/) mt=p
      if (p ~ /\/pt[._]/) pt=p
    }
    /},[[:space:]]*{/ || /\}[[:space:]]*\][[:space:]]*\}/ {
      if (start && species!="") {
        printf "%s\t%s\n", species, pt >> "'"$PT_TSV"'"
        printf "%s\t%s\n", species, mt >> "'"$MT_TSV"'"
      }
      start=0
    }
  ' "$MANIFEST"
fi

# TSV → CSV
awk -F'\t' 'BEGIN{OFS=","} NR==1{print "species","png"; next} {print $1,$2}' "$PT_TSV" >"$PT_CSV"
awk -F'\t' 'BEGIN{OFS=","} NR==1{print "species","png"; next} {print $1,$2}' "$MT_TSV" >"$MT_CSV"

echo "[INFO] First few PT rows:"
sed -n '1,6p' "$PT_CSV"
echo "[INFO] First few MT rows:"
sed -n '1,6p' "$MT_CSV"

# --------- render in R (two pages in one device) ----------
Rscript "${_POLAPLIB_DIR}/scripts/make_png_grid_twopage.R" \
	--pt-list "$PT_CSV" \
	--mt-list "$MT_CSV" \
	--out "$OUTPDF" \
	--pt-title "$PT_TITLE" \
	--mt-title "$MT_TITLE" \
	--rows "$ROWS" \
	--cols "$COLS" \
	--page-width-in "$PW" \
	--page-height-in "$PH" \
	--margin-in "$MG" \
	--gap-mm "$GAP" \
	--label-cex "$LABEL_CEX" \
	--title-cex "$TITLE_CEX" \
	--label-strip-frac "$LABEL_STRIP"

echo "[OK] Wrote: $OUTPDF"
