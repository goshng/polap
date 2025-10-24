#!/usr/bin/env bash
# polap-bash-table-pt.sh
# Version : v0.3.0  (2025-10-15)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Manifest -> TSV + Markdown table summarizing plastid (PT) assemblies.
#
# TSV columns (tab-separated, no quotes, no spaces in header):
#   code  sra  input_gb  n  mean_len  ncbi_acc  ncbi_kb  asm_kb
#
# Markdown columns (8, in this order):
#   Code | SRA | Input Gb | N | Mean len | NCBI (acc) | NCBI kb | Assembly kb
#
set -euo pipefail
IFS=$'\n\t'

MANIFEST="" TSV_OUT="" MD_OUT=""
SORT="code-asc"

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --manifest md/manifest.json --tsv md/pt-summary.tsv --md md/pt-summary.md [--sort code-asc|code-desc|none]

Notes:
  * Two-letter species code is read from .items[].code2 in the manifest.
  * --sort sorts the final TSV rows by 'code' (first column).
EOF
}

have() { command -v "$1" >/dev/null 2>&1; }

# -------------------------- CLI --------------------------
while (($#)); do
	case "$1" in
	--manifest)
		MANIFEST="${2:?}"
		shift 2
		;;
	--tsv)
		TSV_OUT="${2:?}"
		shift 2
		;;
	--md)
		MD_OUT="${2:?}"
		shift 2
		;;
	--sort)
		SORT="${2:?}"
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

[[ -s "${MANIFEST:-}" ]] || {
	echo "[ERR] manifest not found: ${MANIFEST:-<unset>}" >&2
	exit 2
}
[[ -n "${TSV_OUT:-}" ]] || {
	echo "[ERR] --tsv required" >&2
	exit 2
}
[[ -n "${MD_OUT:-}" ]] || {
	echo "[ERR] --md required" >&2
	exit 2
}

mkdir -p "$(dirname "$TSV_OUT")" "$(dirname "$MD_OUT")"

# -------------------------- TSV HEADER --------------------------
# code  sra  input_gb  n  mean_len  ncbi_acc  ncbi_kb  asm_kb
printf "code\tsra\tinput_gb\tn\tmean_len\tncbi_acc\tncbi_kb\tasm_kb\n" >"$TSV_OUT"

# -------------------------- TSV ROWS FROM MANIFEST --------------------------
if ! have jq; then
	echo "[ERR] jq required" >&2
	exit 2
fi

# Emit rows as TSV (code, sra, input_gb, n, mean_len, ncbi_acc, ncbi_kb, asm_kb)
# input_gb prefers data_file_bytes/1e9 else total_bases/1e9
tmp_tsv="$(mktemp -t ptrows.XXXXXX.tsv)"
trap 'rm -f "$tmp_tsv"' EXIT

jq -r '
  .items[] |
  {
    code:   (.code2 // ""),
    sra:    (.data.sra_id // "NA"),
    bytes:  (.data.data_file_bytes // "NA"),
    sum:    (.data.total_bases // "NA"),
    n:      (.data.read_count // "NA"),
    mean:   (.data.mean_length // "NA"),
    acc:    (.pt.ref.ncbi_accession // "NA"),
    nlen:   (.pt.ref.ncbi_len_bp    // "NA"),
    asmlen: (.pt.stats.total_len    // "NA")
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
    .mean,
    .acc,
    ( if (.nlen|tostring)!="NA" then ((.nlen|tonumber)/1e3) else "NA" end ),
    ( if (.asmlen|tostring)!="NA" then ((.asmlen|tonumber)/1e3) else "NA" end )
  ] | @tsv
' "$MANIFEST" >"$tmp_tsv"

# -------------------------- OPTIONAL SORT --------------------------
# Sort by code (first column) asc/desc, keeping header in place
case "$SORT" in
code-asc)
	{
		head -n1 "$TSV_OUT"
		sort -t $'\t' -k1,1 -V "$tmp_tsv"
	} >"${TSV_OUT}.sorted"
	mv "${TSV_OUT}.sorted" "$TSV_OUT"
	;;
code-desc)
	{
		head -n1 "$TSV_OUT"
		sort -t $'\t' -k1,1r -V "$tmp_tsv"
	} >"${TSV_OUT}.sorted"
	mv "${TSV_OUT}.sorted" "$TSV_OUT"
	;;
none | *)
	cat "$tmp_tsv" >>"$TSV_OUT"
	;;
esac

# -------------------------- MARKDOWN TABLE --------------------------
# Format numbers for Markdown, write pretty headers
{
	echo '| Code | SRA | Input Gb | N | Mean len | NCBI (acc) | NCBI kb | Assembly kb |'
	echo '|:--:|:--:|--:|--:|--:|:--|--:|--:|'
} >"$MD_OUT"

# Stream TSV rows (skip header), pretty-print numerics for Markdown
tail -n +2 "$TSV_OUT" | while IFS=$'\t' read -r code sra input_gb n mean_len acc ncbi_kb asm_kb; do
	# pretty numbers
	fmt_gb="$input_gb"
	[[ "$input_gb" != "NA" ]] && fmt_gb=$(awk -v v="$input_gb" 'BEGIN{printf "%.2f", v}')
	fmt_mean="$mean_len"
	[[ "$mean_len" != "NA" ]] && fmt_mean=$(awk -v v="$mean_len" 'BEGIN{printf "%.1f", v}')
	fmt_nkb="$ncbi_kb"
	[[ "$ncbi_kb" != "NA" ]] && fmt_nkb=$(awk -v v="$ncbi_kb" 'BEGIN{printf "%.0f", v}')
	fmt_akb="$asm_kb"
	[[ "$asm_kb" != "NA" ]] && fmt_akb=$(awk -v v="$asm_kb" 'BEGIN{printf "%.0f", v}')

	echo "| ${code} | ${sra} | ${fmt_gb} | ${n} | ${fmt_mean} | ${acc} | ${fmt_nkb} | ${fmt_akb} |" >>"$MD_OUT"
done

echo "[OK] Wrote TSV    : $TSV_OUT"
echo "[OK] Wrote Markdown: $MD_OUT"
