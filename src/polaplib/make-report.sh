#!/usr/bin/env bash
set -euo pipefail

# usage: bash make-report.sh polap-organelle-report.html > report.html
TEMPLATE="${1:-polap-organelle-report.html}"

# helpers --------------------------------------------------------------
esc() { # escape sed-sensitive chars
	printf '%s' "$1" | sed -e 's/[&/\]/\\&/g'
}
read_first() { [ -s "$1" ] && head -n1 "$1" || printf 'NA'; }
sum_fasta_len() {
	# Prefer seqkit (fast); fallback to awk
	local fa="$1"
	if command -v seqkit >/dev/null 2>&1; then
		seqkit stats -Ta "$fa" 2>/dev/null | awk 'NR==2{print $5}' || printf 'NA'
	else
		awk '
      /^>/ {next} {gsub(/[ \t\r\n]/,""); n+=length($0)} END{print n+0}
    ' "$fa" 2>/dev/null || printf 'NA'
	fi
}
human_bp() { # like 122725 -> 122.7 kbp ; 100751943 -> 100.8 Mbp
	awk -v x="$1" 'BEGIN{
    if (x=="" || x==0) {print "NA"; exit}
    if (x>=1e9)  printf "%.1f Gbp", x/1e9;
    else if (x>=1e6) printf "%.1f Mbp", x/1e6;
    else if (x>=1e3) printf "%.1f kbp", x/1e3;
    else printf "%d bp", x+0;
  }'
}
# values ---------------------------------------------------------------
RUN_ROOT="."
RUN_NAME="$(basename "$RUN_ROOT")"
RUN_DATE="$(date +%F)"

POLAP_VERSION="$(grep -m1 'POLAP:' "$RUN_ROOT/polap.log" 2>/dev/null | awk '{print $2}' || true)"
POLAP_VERSION="${POLAP_VERSION:-NA}"

POLAP_CMDLINE="$(grep -m1 '^\\[.*\\] CMD:' "$RUN_ROOT/polap.log" 2>/dev/null | sed 's/.*CMD:[ ]*//; s/[[:space:]]*$//' || true)"
POLAP_CMDLINE="${POLAP_CMDLINE:-NA}"

READ_BASES_TOTAL="$(read_first "$RUN_ROOT/annotate-read-pt/pt1/01-contig/l.txt")"
READ_BASES_SUBSAMPLED="$(read_first "$RUN_ROOT/annotate-read-pt/pt1/01-contig/l.subsample.txt")"

# PT assembled size (sum of graph_final.fasta in pt1 if present, else pt0)
PT_FASTA=""
if [ -s "$RUN_ROOT/annotate-read-pt/pt1/30-contigger/graph_final.fasta" ]; then
	PT_FASTA="$RUN_ROOT/annotate-read-pt/pt1/30-contigger/graph_final.fasta"
elif [ -s "$RUN_ROOT/annotate-read-pt/pt/30-contigger/graph_final.fasta" ]; then
	PT_FASTA="$RUN_ROOT/annotate-read-pt/pt/30-contigger/graph_final.fasta"
fi
PT_SIZE="NA"
if [ -n "$PT_FASTA" ]; then
	PT_SUM="$(sum_fasta_len "$PT_FASTA")"
	PT_SIZE="$(human_bp "$PT_SUM")"
fi

# MT rounds present (mt0, mt1, mt2) by presence of graph_final.gfa
MT_ROUNDS=""
for r in mt0 mt1 mt2; do
	if [ -s "$RUN_ROOT/mtseed/$r/30-contigger/graph_final.gfa" ]; then
		MT_ROUNDS="$MT_ROUNDS $r"
	fi
done
MT_ROUNDS="$(echo "$MT_ROUNDS" | xargs || printf 'none')"

# placeholder replacement ---------------------------------------------
sed -e "s|%%RUN_NAME%%|$(esc "$RUN_NAME")|g" \
	-e "s|%%RUN_DATE%%|$(esc "$RUN_DATE")|g" \
	-e "s|%%POLAP_VERSION%%|$(esc "$POLAP_VERSION")|g" \
	-e "s|%%POLAP_CMDLINE%%|$(esc "$POLAP_CMDLINE")|g" \
	-e "s|%%READ_BASES_TOTAL%%|$(esc "$(human_bp "$READ_BASES_TOTAL")")|g" \
	-e "s|%%READ_BASES_SUBSAMPLED%%|$(esc "$(human_bp "$READ_BASES_SUBSAMPLED")")|g" \
	-e "s|%%PT_SIZE%%|$(esc "$PT_SIZE")|g" \
	-e "s|%%MT_ROUNDS%%|$(esc "$MT_ROUNDS")|g" \
	"$TEMPLATE"
