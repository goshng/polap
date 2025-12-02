#!/usr/bin/env bash
# polap-bash-hifi-ptmt-sheet.sh
# Version : v1.1.0  (2025-12-02)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Draw mitogenome assembly graphs from 4 pipelines (pmat, tippo, himt, oatk)
# per species, per row, using v5-0-auto-manifest.json.
#
# For each species, columns are:
#   1. PMAT (pmat2 in JSON)
#   2. TIPPo
#   3. HiMT
#   4. Oatk (prefers c30, then c20, then c10, then mt_png)
#
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
MANIFEST=""
OUTPDF=""
PW=11 PH=8.5 MG=0.45 GAP=4.0
LABEL_CEX=0.55 TITLE_CEX=1.0
PTMT_TITLE="Mitogenome assembly graphs (PMAT, TIPPo, HiMT, Oatk)"
SPECIES_CODES_TXT="${_POLAPLIB_DIR}/polap-species-codes-hifi.txt"

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --manifest v5-0-auto-manifest.json --out md/hifi-mt-4pipelines.pdf \\
                   --species-codes polap-species-codes-hifi.txt \\
                   [--base-dir <dir>] \\
                   [--page-width-in IN] [--page-height-in IN] [--margin-in IN] \\
                   [--gap-mm MM] [--label-cex CEX] [--title-cex CEX] \\
                   [--title "Mitogenome assembly graphs (PMAT, TIPPo, HiMT, Oatk)"]

Required:
  --manifest        v5-0-auto-manifest.json for this analysis
  --out             Output PDF path

Pipelines (columns, in order):
  1. pmat  (JSON key: pmat2, field: mt_png)
  2. tippo (JSON key: tippo, field: mt_png)
  3. himt  (JSON key: himt,  field: mt_png)
  4. oatk  (JSON key: oatk,  fields: c30_png, c20_png, c10_png, mt_png)
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
	--title)
		PTMT_TITLE="${2:?}"
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
	# BASE_DIR is typically something like Species/v5/0/...
	# We go one level up so relative PNG paths still work as in v5-0-auto-manifest.json.
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
# Build unordered list from v5-0-auto-manifest.json
#   columns: species, pmat_png, tippo_png, himt_png, oatk_png
# -----------------------------------------------------------------------------#
echo -e "species\tpmat_png\ttippo_png\thimt_png\toatk_png" >"$TSV_UN"

if command -v jq >/dev/null 2>&1; then
	jq -r --arg SR "$SR" '
	  .items[]
	  | . as $it
	  | $it.species as $sp
	  | ($it.pmat2.mt_png // "") as $pmat_png
	  | ($it.tippo.mt_png // "") as $tippo_png
	  | ($it.himt.mt_png  // "") as $himt_png
	  | (
	      $it.oatk.c30_png
	      // $it.oatk.c20_png
	      // $it.oatk.c10_png
	      // $it.oatk.mt_png
	      // ""
	    ) as $oatk_png
	  | [
	      $sp,
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

# -----------------------------------------------------------------------------#
# Reorder to match species-codes.txt
# -----------------------------------------------------------------------------#
reorder_tsv() {
	local order="$1" unordered="$2" out="$3"
	echo -e "species\tpmat_png\ttippo_png\thimt_png\toatk_png" >"$out"
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
awk -F'\t' 'BEGIN{OFS=","} NR==1{print "species","pmat_png","tippo_png","himt_png","oatk_png"; next} {print $1,$2,$3,$4,$5}' \
	"$TSV" >"$CSV"

# -----------------------------------------------------------------------------#
# Render with R
# -----------------------------------------------------------------------------#
Rscript "${_POLAPLIB_DIR}/scripts/polap-r-hifi-ptmt-sheet.R" \
	--list "$CSV" \
	--out "$OUTPDF" \
	--title "$PTMT_TITLE" \
	--page-width-in "$PW" --page-height-in "$PH" --margin-in "$MG" \
	--gap-mm "$GAP" --label-cex "$LABEL_CEX" --title-cex "$TITLE_CEX"

echo "[OK] Wrote: $OUTPDF"
