#!/usr/bin/env bash
# polap-bash-table-mt.sh
# Version : v0.3.0  (2025-10-15)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Mitochondrial genome summary: manifest -> TSV -> Markdown table
# Output TSV uses tabs (no commas, no quotes).
#
# TSV columns:
#   code	sra	input_gb	mean_read_len_bases	ncbi_acc	ncbi_len_bp	asm_segments	asm_bases
#
# Markdown columns:
#   Code | SRA | Input Gb | Mean read length (bases) | NCBI accession | NCBI length | Assembly segments | Assembly bases
#
set -euo pipefail
IFS=$'\n\t'

MANIFEST="" TSV_OUT="" MD_OUT=""

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --manifest md/manifest.json --tsv md/mt-summary.tsv --md md/mt-summary.md
EOF
}

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
[[ -n "${TSV_OUT:-}" ]] || {
	echo "[ERR] --tsv required" >&2
	exit 2
}
[[ -n "${MD_OUT:-}" ]] || {
	echo "[ERR] --md required" >&2
	exit 2
}

mkdir -p "$(dirname "$TSV_OUT")" "$(dirname "$MD_OUT")"

have() { command -v "$1" >/dev/null 2>&1; }

# ------------------------------------------------------------------------------
# TSV (tab-separated, no quotes)
# ------------------------------------------------------------------------------
{
	echo -e "code\tsra\tinput_gb\tmean_read_len_bases\tncbi_acc\tncbi_len_bp\tasm_segments\tasm_bases"
} >"$TSV_OUT"

if have jq; then
	jq -r '
    .items[] |
    {
      code:   (.code2 // ""),
      sra:    (.data.sra_id // "NA"),
      bytes:  (.data.data_file_bytes // "NA"),
      sum:    (.data.total_bases // "NA"),
      meanlen:(.data.mean_length // "NA"),
      acc:    (.mt.ref.ncbi_accession // "NA"),
      nlen:   (.mt.ref.ncbi_len_bp // "NA"),
      nseg:   (.mt.stats.n_segments // "NA"),
      asmlen: (.mt.stats.total_len // "NA")
    } |
    [
      .code,
      .sra,
      ( if (.bytes|tostring)!="NA" and (.bytes|tonumber)>0
        then ((.bytes|tonumber)/1e9)
        else ( if (.sum|tostring)!="NA" and (.sum|tonumber)>0
               then ((.sum|tonumber)/1e9) else "NA" end )
        end ),
      .meanlen,
      .acc,
      .nlen,
      .nseg,
      .asmlen
    ] | @tsv
  ' "$MANIFEST" >>"$TSV_OUT"
else
	echo "[ERR] jq required" >&2
	exit 2
fi

# ------------------------------------------------------------------------------
# Markdown (human-readable)
# ------------------------------------------------------------------------------
{
	echo '| Code | SRA | Input Gb | Mean read length (bases) | NCBI accession | NCBI length | Assembly segments | Assembly bases |'
	echo '|:--:|:--:|--:|--:|:--|--:|--:|--:|'
} >"$MD_OUT"

tail -n +2 "$TSV_OUT" | while IFS=$'\t' read -r code sra input_gb meanlen acc nlen nseg asmbases; do
	# numeric formatting
	fmt_gb="$input_gb"
	[[ "$input_gb" != "NA" ]] && fmt_gb=$(awk -v v="$input_gb" 'BEGIN{printf "%.2f", v}')
	fmt_mean="$meanlen"
	[[ "$meanlen" != "NA" ]] && fmt_mean=$(awk -v v="$meanlen" 'BEGIN{printf "%.1f", v}')
	fmt_nlen="$nlen"
	[[ "$nlen" != "NA" ]] && fmt_nlen=$(awk -v v="$nlen" 'BEGIN{printf "%.0f", v}')
	fmt_nseg="$nseg"
	[[ "$nseg" != "NA" ]] && fmt_nseg=$(awk -v v="$nseg" 'BEGIN{printf "%.0f", v}')
	fmt_asmb="$asmbases"
	[[ "$asmbases" != "NA" ]] && fmt_asmb=$(awk -v v="$asmbases" 'BEGIN{printf "%.0f", v}')

	echo "| ${code} | ${sra} | ${fmt_gb} | ${fmt_mean} | ${acc} | ${fmt_nlen} | ${fmt_nseg} | ${fmt_asmb} |" >>"$MD_OUT"
done

echo "[OK] Wrote TSV    : $TSV_OUT"
echo "[OK] Wrote Markdown: $MD_OUT"
