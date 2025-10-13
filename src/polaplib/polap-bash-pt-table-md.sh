#!/usr/bin/env bash
# polap-bash-pt-table-md.sh
# Version : v0.2.1  (2025-10-14)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Manifest → Markdown table summarizing plastid assemblies.
# Columns (exactly 8, in this order):
#   Code | SRA | Input Gb | N | Mean len | NCBI (acc) | NCBI kb | Assembly kb
#
set -euo pipefail
IFS=$'\n\t'

MANIFEST="" OUT=""
usage() {
	cat <<EOF
Usage:
  $(basename "$0") --manifest md/manifest.json --out md/pt-summary.md
  (Two-letter species code is read from .items[].code2 in the manifest.)
EOF
}

while (($#)); do
	case "$1" in
	--manifest)
		MANIFEST="${2:?}"
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
		usage
		exit 2
		;;
	esac
done

[[ -s "$MANIFEST" ]] || {
	echo "[ERR] manifest not found: $MANIFEST" >&2
	exit 2
}
mkdir -p "$(dirname "$OUT")"

have() { command -v "$1" >/dev/null 2>&1; }

# header (8 columns)
{
	echo '| Code | SRA | Input Gb | N | Mean len | NCBI (acc) | NCBI kb | Assembly kb |'
	echo '|:--:|:--:|--:|--:|--:|:--|--:|--:|'
} >"$OUT"

if have jq; then
	# Emit a TSV with exactly these 8 fields (in order):
	# code, sra, gb, n, al, acc, kb, asm
	jq -r '
	  .items[] |
	  {
	    code: (.code2 // ""),
	    sra:  (.data.sra_id // "NA"),
	    bytes:(.data.data_file_bytes // "NA"),
	    sum:  (.data.total_bases // "NA"),
	    n:    (.data.read_count // "NA"),
	    al:   (.data.mean_length // "NA"),
	    acc:  (.pt.ref.ncbi_accession // "NA"),
	    nlen: (.pt.ref.ncbi_len_bp // "NA"),
	    asml: (.pt.stats.total_len // "NA")
	  } |
	  [
	    .code,
	    .sra,
	    ( if (.bytes|tostring)!="NA" and (.bytes|tonumber)>0
	      then ((.bytes|tonumber)/1e9)
	      else ( if (.sum|tostring)!="NA" and (.sum|tonumber)>0
	             then ((.sum|tonumber)/1e9) else "NA" end )
	      end ),
	    .n,
	    .al,
	    .acc,
	    ( if (.nlen|tostring)!="NA" then ((.nlen|tonumber)/1e3) else "NA" end ),
	    ( if (.asml|tostring)!="NA" then ((.asml|tonumber)/1e3) else "NA" end )
	  ] | @tsv
	' "$MANIFEST" |
		while IFS=$'\t' read -r code sra gb n al acc kb asm; do
			# pretty numbers
			fmt_gb="$gb"
			[[ "$gb" != "NA" ]] && fmt_gb=$(awk -v v="$gb" 'BEGIN{printf "%.2f", v}')
			fmt_al="$al"
			[[ "$al" != "NA" ]] && fmt_al=$(awk -v v="$al" 'BEGIN{printf "%.1f", v}')
			fmt_kb="$kb"
			[[ "$kb" != "NA" ]] && fmt_kb=$(awk -v v="$kb" 'BEGIN{printf "%.0f", v}')
			fmt_asm="$asm"
			[[ "$asm" != "NA" ]] && fmt_asm=$(awk -v v="$asm" 'BEGIN{printf "%.0f", v}')
			echo "| ${code} | ${sra} | ${fmt_gb} | ${n} | ${fmt_al} | ${acc} | ${fmt_kb} | ${fmt_asm} |" >>"$OUT"
		done
else
	echo "[ERR] jq required" >&2
	exit 2
fi

echo "[OK] Wrote Markdown: $OUT"
