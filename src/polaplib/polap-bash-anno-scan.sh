#!/usr/bin/env bash
# polap-bash-anno-scan.sh
# Version : v1.1.2  (2025-10-07)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# From a manifest JSON, produce:
#   <out-dir>/anno-mt.csv  (species,mt_genes,pt_genes,genome_kb)
#   <out-dir>/anno-pt.csv  (species,mt_genes,pt_genes,genome_kb)
#
# Reads contig-annotation-depth-table.txt files, summing Length, MT, PT columns.

set -euo pipefail
IFS=$'\n\t'

AWK_CORE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/scripts/polap-awk-anno-sum.awk"
[[ -f "$AWK_CORE" ]] || {
	echo "[ERR] Missing AWK core: $AWK_CORE" >&2
	exit 2
}

BASE_DIR=""
MANIFEST="" OUT_DIR="." TIER="" INUM=""
usage() {
	cat <<EOF
Usage:
  $(basename "$0") --manifest md/manifest.json --out-dir md [--tier v6] [--inum 0]
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
	--out-dir)
		OUT_DIR="${2:?}"
		shift 2
		;;
	--tier)
		TIER="${2:?}"
		shift 2
		;;
	--inum)
		INUM="${2:?}"
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
mkdir -p "$OUT_DIR"

SR="$(cd "$(dirname "$BASE_DIR")" && pwd)"

MT_OUT="$OUT_DIR/anno-mt.csv"
PT_OUT="$OUT_DIR/anno-pt.csv"
echo "species,mt_genes,pt_genes,genome_kb" >"$MT_OUT"
echo "species,mt_genes,pt_genes,genome_kb" >"$PT_OUT"

have() { command -v "$1" >/dev/null 2>&1; }

# fallback: total length from a GFA
gfa_total_len() {
	local gfa="$1"
	[[ -s "$gfa" ]] || {
		echo "NA"
		return 0
	}
	awk '
    $1=="S"{
      ln=0
      if($3!="*") ln=length($3)
      for(i=4;i<=NF;i++) if($i~/^LN:i:/){ split($i,a,":"); ln=a[3]; break }
      sum+=ln
    } END{ print (sum?sum:"NA") }
  ' "$gfa"
}

# species list
if have jq; then
	mapfile -t SP < <(jq -r '.items[].species' "$MANIFEST" | sort -u)
else
	mapfile -t SP < <(awk -F'"' '/"species"[[:space:]]*:/ {print $4}' "$MANIFEST" | sort -u)
fi

# optional lengths from manifest
declare -A LEN_MT LEN_PT
if have jq; then
	while IFS=$'\t' read -r s l; do [[ -n "$s" ]] && LEN_MT["$s"]="$l"; done \
		< <(jq -r '.items[] | select(.mt.stats.total_len?) | [.species, .mt.stats.total_len] | @tsv' "$MANIFEST")
	while IFS=$'\t' read -r s l; do [[ -n "$s" ]] && LEN_PT["$s"]="$l"; done \
		< <(jq -r '.items[] | select(.pt.stats.total_len?) | [.species, .pt.stats.total_len] | @tsv' "$MANIFEST")
fi

for sp in "${SP[@]}"; do
	base="$SR/$sp"
	tier="${TIER:-$(jq -r --arg s "$sp" '.items[] | select(.species==$s) | .tier' "$MANIFEST" 2>/dev/null | head -n1)}"
	inum="${INUM:-$(jq -r --arg s "$sp" '.items[] | select(.species==$s) | .inum' "$MANIFEST" 2>/dev/null | head -n1)}"
	[[ -z "$tier" || "$tier" == "null" ]] && tier="v6"
	[[ -z "$inum" || "$inum" == "null" ]] && inum="0"
	run="$base/$tier/$inum/polap-readassemble"

	########################################################################
	# MT table
	########################################################################
	mt_tab="$run/mtseed/mt1/contig-annotation-depth-table.txt"
	[[ ! -s "$mt_tab" ]] && mt_tab="$run/contig-annotation-depth-table.txt"

	mt_genes="NA"
	pt_in_mt="NA"
	mt_kb="NA"

	if [[ -s "$mt_tab" ]]; then
		mt_out=$(awk -f "$AWK_CORE" "$mt_tab" 2>/dev/null || true)
		if [[ -n "${mt_out:-}" ]]; then
			IFS=$'\t' read -r len_bp len_kb mt_sum pt_sum nrow <<<"$mt_out"
			len_bp="${len_bp:-0}"
			mt_sum="${mt_sum:-0}"
			pt_sum="${pt_sum:-0}"
			mt_genes="$mt_sum"
			pt_in_mt="$pt_sum"
			mt_kb="$(awk -v v="$len_bp" 'BEGIN{printf "%.0f", v/1000.0}')"
		else
			echo "[WARN] Empty or unparsable MT table: $mt_tab" >&2
		fi
	else
		len=$(gfa_total_len "$run/mt.1.gfa")
		[[ "$len" == "NA" ]] && len=$(gfa_total_len "$run/mt1/assembly_graph.gfa")
		[[ "$len" != "NA" ]] && mt_kb="$(awk -v v="$len" 'BEGIN{printf "%.0f", v/1000.0}')"
	fi
	echo "$sp,$mt_genes,$pt_in_mt,$mt_kb" >>"$MT_OUT"

	########################################################################
	# PT table
	########################################################################
	pt_tab="$run/annotate-read-pt/pt1/pt-contig-annotation-depth-table.txt"
	[[ ! -s "$pt_tab" ]] && pt_tab="$run/pt/pt-contig-annotation-depth-table.txt"

	pt_genes="NA"
	mt_in_pt="NA"
	pt_kb="NA"

	if [[ -s "$pt_tab" ]]; then
		pt_out=$(awk -f "$AWK_CORE" "$pt_tab" 2>/dev/null || true)
		if [[ -n "${pt_out:-}" ]]; then
			IFS=$'\t' read -r len_bp len_kb mt_sum pt_sum nrow <<<"$pt_out"
			len_bp="${len_bp:-0}"
			mt_sum="${mt_sum:-0}"
			pt_sum="${pt_sum:-0}"
			mt_in_pt="$mt_sum"
			pt_genes="$pt_sum"
			pt_kb="$(awk -v v="$len_bp" 'BEGIN{printf "%.0f", v/1000.0}')"
		else
			echo "[WARN] Empty or unparsable PT table: $pt_tab" >&2
		fi
	else
		len=$(gfa_total_len "$run/pt.1.gfa")
		[[ "$len" == "NA" ]] && len=$(gfa_total_len "$run/pt-pt.1.gfa")
		[[ "$len" == "NA" ]] && len=$(gfa_total_len "$run/pt1/assembly_graph.gfa")
		[[ "$len" != "NA" ]] && pt_kb="$(awk -v v="$len" 'BEGIN{printf "%.0f", v/1000.0}')"
	fi
	echo "$sp,$mt_in_pt,$pt_genes,$pt_kb" >>"$PT_OUT"
done

echo "[OK] MT annotations: $MT_OUT"
echo "[OK] PT annotations: $PT_OUT"
