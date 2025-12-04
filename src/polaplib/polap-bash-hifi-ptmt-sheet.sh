#!/usr/bin/env bash
# polap-bash-hifi-ptmt-sheet.sh
# Version : v1.1.3  (2025-12-02)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Draw mitogenome assembly graphs from 4 pipelines (pmat, tippo, himt, oatk)
# per species, per row, using v5-0-auto-manifest.json.
#
set -euo pipefail
IFS=$'\n\t'

_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export _POLAPLIB_DIR

# -----------------------------------------------------------------------------#
# Defaults and CLI args
# -----------------------------------------------------------------------------#
BASE_DIR=""
MANIFEST=""
OUTPDF=""
PW=11 PH=8.5 MG=0.45 GAP=4.0
LABEL_CEX=0.55 TITLE_CEX=1.0
PTMT_TITLE="Mitogenome assembly graphs (PMAT, TIPPo, HiMT, Oatk)"
SUBTITLE=""
ROWS_PER_PAGE=0 # 0 = auto (= all on one page)
START_PAGE=1
SPECIES_CODES_TXT="${_POLAPLIB_DIR}/polap-species-codes-hifi.txt"

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --manifest v5-0-auto-manifest.json --out md/hifi-mt-4pipelines.pdf \\
                   --species-codes polap-species-codes-hifi.txt \\
                   [--base-dir <dir>] \\
                   [--page-width-in IN] [--page-height-in IN] [--margin-in IN] \\
                   [--gap-mm MM] [--label-cex CEX] [--title-cex CEX] \\
                   [--rows-per-page N] [--start-page N] \\
                   [--title "text"] [--subtitle "text"]

Required:
  --manifest        v5-0-auto-manifest.json for this analysis
  --out             Output PDF path
EOF
}

while (($#)); do
	case "$1" in
	--species-codes)
		SPECIES_CODES_TXT="${2:?}"
		shift 2
		;;
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
	--rows-per-page)
		ROWS_PER_PAGE="${2:?}"
		shift 2
		;;
	--start-page)
		START_PAGE="${2:?}"
		shift 2
		;;
	--title)
		PTMT_TITLE="${2:?}"
		shift 2
		;;
	--subtitle)
		SUBTITLE="${2:?}"
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
TSV_UN="$(mktemp -t hifi-mt4pip.unordered.XXXXXX.tsv)"
TSV="$(mktemp -t hifi-mt4pip.ordered.XXXXXX.tsv)"
CSV="$(mktemp -t hifi-mt4pip.XXXXXX.csv)"
ORDER_FILE="$(mktemp -t hifi-species.order.XXXXXX.txt)"
trap 'rm -f "$TSV_UN" "$TSV" "$CSV" "$ORDER_FILE"' EXIT

# -----------------------------------------------------------------------------#
# Species order file from species-codes.txt
# -----------------------------------------------------------------------------#
if [[ -s "$SPECIES_CODES_TXT" ]]; then
	awk 'NR==1{next} NF>=2 {print $2}' "$SPECIES_CODES_TXT" >"$ORDER_FILE"
fi

# -----------------------------------------------------------------------------#
# Build unordered list from manifest
#   columns: species, code2, pmat_png, tippo_png, himt_png, oatk_png
# -----------------------------------------------------------------------------#
echo -e "species\tcode2\tpmat_png\ttippo_png\thimt_png\toatk_png" >"$TSV_UN"

if command -v jq >/dev/null 2>&1; then
	jq -r --arg SR "$SR" '
	  .items[]
	  | . as $it
	  | $it.species as $sp
	  | $it.code2   as $code2
	  | ($it.pmat2.mt_png // "") as $pmat_png
	  | ($it.tippo.mt_png // "") as $tippo_png
	  | ($it.himt.mt_png  // "") as $himt_png
	  | ($it.oatk.c30_png  // "") as $oatk_png
	  | [
	      $sp,
	      $code2,
	      (if $pmat_png=="" then "" else ($SR + "/" + $pmat_png) end),
	      (if $tippo_png=="" then "" else ($SR + "/" + $tippo_png) end),
	      (if $himt_png==""  then "" else ($SR + "/" + $himt_png)  end),
	      (if $oatk_png==""  then "" else ($SR + "/" + $oatk_png)  end)
	    ]
	  | @tsv
	' "$MANIFEST" >>"$TSV_UN"
else
	echo "[ERR] jq required" >&2
	exit 2
fi

# if command -v jq >/dev/null 2>&1; then
# 	jq -r --arg SR "$SR" '
# 	  .items[]
# 	  | . as $it
# 	  | $it.species as $sp
# 	  | $it.code2   as $code2
# 	  | ($it.pmat2.mt_png // "") as $pmat_png
# 	  | ($it.tippo.mt_png // "") as $tippo_png
# 	  | ($it.himt.mt_png  // "") as $himt_png
# 	  | (
# 	      $it.oatk.c30_png
# 	      // $it.oatk.c20_png
# 	      // $it.oatk.c10_png
# 	      // $it.oatk.mt_png
# 	      // ""
# 	    ) as $oatk_png
# 	  | [
# 	      $sp,
# 	      $code2,
# 	      (if $pmat_png=="" then "" else ($SR + "/" + $pmat_png) end),
# 	      (if $tippo_png=="" then "" else ($SR + "/" + $tippo_png) end),
# 	      (if $himt_png==""  then "" else ($SR + "/" + $himt_png)  end),
# 	      (if $oatk_png==""  then "" else ($SR + "/" + $oatk_png)  end)
# 	    ]
# 	  | @tsv
# 	' "$MANIFEST" >>"$TSV_UN"
# else
# 	echo "[ERR] jq required" >&2
# 	exit 2
# fi

# -----------------------------------------------------------------------------#
# Reorder to match species-codes.txt (key = species)
# -----------------------------------------------------------------------------#
reorder_tsv() {
	local order="$1" unordered="$2" out="$3"
	echo -e "species\tcode2\tpmat_png\ttippo_png\thimt_png\toatk_png" >"$out"
	if [[ -s "$order" ]]; then
		awk -F'\t' -v OFS='\t' '
		  NR==FNR {
		    if($0!~/^[[:space:]]*$/){
		      n++; ord[n]=$0; want[$0]=1
		    }
		    next
		  }
		  FNR==1 { next }
		  { lines[$1]=$0 }
		  END {
		    for(i=1;i<=n;i++){
		      s=ord[i]
		      if(s in lines) print lines[s]
		    }
		    for(s in lines) if(!(s in want)) print lines[s]
		  }
		' "$order" "$unordered" >>"$out"
	else
		tail -n +2 "$unordered" >>"$out"
	fi
}

reorder_tsv "$ORDER_FILE" "$TSV_UN" "$TSV"

# -----------------------------------------------------------------------------#
# Convert to CSV for the R script
# -----------------------------------------------------------------------------#
awk -F'\t' 'BEGIN{OFS=","} NR==1{print "species","code2","pmat_png","tippo_png","himt_png","oatk_png"; next} {print $1,$2,$3,$4,$5,$6}' \
	"$TSV" >"$CSV"

# -----------------------------------------------------------------------------#
# Render with R
# -----------------------------------------------------------------------------#
Rscript "${_POLAPLIB_DIR}/scripts/polap-r-hifi-ptmt-sheet.R" \
	--list "$CSV" \
	--out "$OUTPDF" \
	--title "$PTMT_TITLE" \
	--subtitle "$SUBTITLE" \
	--rows-per-page "$ROWS_PER_PAGE" \
	--start-page "$START_PAGE" \
	--page-width-in "$PW" --page-height-in "$PH" --margin-in "$MG" \
	--gap-mm "$GAP" --label-cex "$LABEL_CEX" --title-cex "$TITLE_CEX"

cp "$CSV" "${OUTPDF%.pdf}.csv"

echo "[OK] Wrote: ${OUTPDF%.pdf}.csv"
echo "[OK] Wrote: $OUTPDF"
