#!/usr/bin/env bash
# polap-bash-ptmt-sheet.sh
# Version : v1.2.0  (2025-10-14)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Build two-page PDF: page 1 = PT, page 2 = MT, from a manifest JSON.
# Ordered by species sequence in species-codes.txt (code species, space-delimited).
# Accepts an explicit --base-dir to set the working root for relative PNG paths.
#
set -euo pipefail
IFS=$'\n\t'

_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export _POLAPLIB_DIR

# -----------------------------------------------------------------------------#
# Defaults and CLI args
# -----------------------------------------------------------------------------#
BASE_DIR=""
MANIFEST="" OUTPDF=""
ANN_PT="" ANN_MT=""
ROWS=0 COLS=0 PW=11 PH=8.5 MG=0.45 GAP=4.0
LABEL_CEX=0.55 TITLE_CEX=1.0 LABEL_STRIP=0.16
PT_TITLE="Plastid assemblies (pt)"
MT_TITLE="Mitochondrial assemblies (mt)"
CODES_FILE="${_POLAPLIB_DIR}/species-codes.txt" # code species (space-delimited)

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --manifest md/manifest.json --out out.pdf
                   --annot-pt md/anno-pt.csv --annot-mt md/anno-mt.csv
                   [--base-dir <dir>]
                   [--rows N] [--cols N]
                   [--page-width-in IN] [--page-height-in IN] [--margin-in IN]
                   [--gap-mm MM] [--label-cex CEX] [--title-cex CEX] [--label-strip-frac FRAC]
                   [--pt-title "text"] [--mt-title "text"]
EOF
}

while (($#)); do
	case "$1" in
	--base-dir)
		BASE_DIR="${2:?}"
		shift 2
		;;
	--manifest)
		MANIFEST="${2:?}"
		shift 2
		;;
	--out)
		OUTPDF="${2:?}"
		shift 2
		;;
	--annot-pt)
		ANN_PT="${2:?}"
		shift 2
		;;
	--annot-mt)
		ANN_MT="${2:?}"
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

[[ -s "$MANIFEST" && -n "$OUTPDF" ]] || {
	echo "[ERR] need --manifest and --out" >&2
	exit 2
}
[[ -s "$ANN_PT" && -s "$ANN_MT" ]] || {
	echo "[ERR] need --annot-pt and --annot-mt CSVs" >&2
	exit 2
}

mkdir -p "$(dirname "$OUTPDF")"

# -----------------------------------------------------------------------------#
# Determine species root (SR)
# -----------------------------------------------------------------------------#
if [[ -n "$BASE_DIR" ]]; then
	SR="$(cd "$(dirname "$BASE_DIR")" && pwd)"
else
	SR="$(cd "$(dirname "$MANIFEST")/.." && pwd)"
fi

# -----------------------------------------------------------------------------#
# Temp files
# -----------------------------------------------------------------------------#
PT_TSV_UN="$(mktemp -t ptlist.unordered.XXXXXX.tsv)"
MT_TSV_UN="$(mktemp -t mtlist.unordered.XXXXXX.tsv)"
PT_TSV="$(mktemp -t ptlist.ordered.XXXXXX.tsv)"
MT_TSV="$(mktemp -t mtlist.ordered.XXXXXX.tsv)"
PT_CSV="$(mktemp -t ptlist.XXXXXX.csv)"
MT_CSV="$(mktemp -t mtlist.XXXXXX.csv)"
ORDER_FILE="$(mktemp -t species.order.XXXXXX.txt)"
trap 'rm -f "$PT_TSV_UN" "$MT_TSV_UN" "$PT_TSV" "$MT_TSV" "$PT_CSV" "$MT_CSV" "$ORDER_FILE"' EXIT

# -----------------------------------------------------------------------------#
# Species order file from species-codes.txt
# -----------------------------------------------------------------------------#
if [[ -s "$CODES_FILE" ]]; then
	awk 'NR==1{next} NF>=2 {print $2}' "$CODES_FILE" >"$ORDER_FILE"
fi

# -----------------------------------------------------------------------------#
# Build unordered lists from manifest
# -----------------------------------------------------------------------------#
echo -e "species\tpng" >"$PT_TSV_UN"
echo -e "species\tpng" >"$MT_TSV_UN"

if command -v jq >/dev/null 2>&1; then
	jq -r --arg SR "$SR" '
    .items[] | {species} + {png: (.pt.png // .pt.file.png // "")} |
    [ .species, (if .png=="" then "" else ($SR + "/" + .png) end) ] | @tsv
  ' "$MANIFEST" >>"$PT_TSV_UN"

	jq -r --arg SR "$SR" '
    .items[] | {species} + {png: (.mt.png // .mt.file.png // "")} |
    [ .species, (if .png=="" then "" else ($SR + "/" + .png) end) ] | @tsv
  ' "$MANIFEST" >>"$MT_TSV_UN"
else
	echo "[ERR] jq required" >&2
	exit 2
fi

# -----------------------------------------------------------------------------#
# Reorder to match species-codes.txt
# -----------------------------------------------------------------------------#
reorder_tsv() {
	local order="$1" unordered="$2" out="$3"
	echo -e "species\tpng" >"$out"
	if [[ -s "$order" ]]; then
		awk -F'\t' -v OFS='\t' '
      NR==FNR { if($0!~/^[[:space:]]*$/){ n++; ord[n]=$0; want[$0]=1 } next }
      FNR==1 { next }
      { lines[$1]=$0 }
      END {
        for(i=1;i<=n;i++){ s=ord[i]; if(s in lines) print lines[s] }
        for(s in lines) if(!(s in want)) print lines[s]
      }
    ' "$order" "$unordered" >>"$out"
	else
		tail -n +2 "$unordered" >>"$out"
	fi
}

reorder_tsv "$ORDER_FILE" "$PT_TSV_UN" "$PT_TSV"
reorder_tsv "$ORDER_FILE" "$MT_TSV_UN" "$MT_TSV"

awk -F'\t' 'BEGIN{OFS=","} NR==1{print "species","png"; next} {print $1,$2}' "$PT_TSV" >"$PT_CSV"
awk -F'\t' 'BEGIN{OFS=","} NR==1{print "species","png"; next} {print $1,$2}' "$MT_TSV" >"$MT_CSV"

# -----------------------------------------------------------------------------#
# Render with R
# -----------------------------------------------------------------------------#
Rscript "${_POLAPLIB_DIR}/scripts/polap-r-ptmt-sheet.R" \
	--pt-list "$PT_CSV" --mt-list "$MT_CSV" \
	--annot-pt "$ANN_PT" --annot-mt "$ANN_MT" \
	--out "$OUTPDF" \
	--pt-title "$PT_TITLE" --mt-title "$MT_TITLE" \
	--rows "$ROWS" --cols "$COLS" \
	--page-width-in "$PW" --page-height-in "$PH" --margin-in "$MG" \
	--gap-mm "$GAP" --label-cex "$LABEL_CEX" --title-cex "$TITLE_CEX" \
	--label-strip-frac "$LABEL_STRIP"

echo "[OK] Wrote: $OUTPDF"
