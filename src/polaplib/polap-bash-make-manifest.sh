#!/usr/bin/env bash
# polaplib/polap-bash-make-manifest.sh
# Version : v0.2.1  (hardened)
set -euo pipefail
IFS=$'\n\t'

# auto base
_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export _POLAPLIB_DIR

# defaults
SET="some"
TIER="v6"
INUM="0"
OUT="md/manifest.json"
INCLUDE_PTPT=0
PRETTY=0
while (($#)); do
	case "$1" in
	--set)
		SET="${2:?}"
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
	--out)
		OUT="${2:?}"
		shift 2
		;;
	--include-ptpt)
		INCLUDE_PTPT=1
		shift
		;;
	--pretty)
		PRETTY=1
		shift
		;;
	*)
		echo "[ERR] Unknown arg: $1" >&2
		exit 2
		;;
	esac
done
mkdir -p "$(dirname "$OUT")"

have() { command -v "$1" >/dev/null 2>&1; }

CSV_CFG="${_POLAPLIB_DIR}/../polap-bash-report.csv" # adjust if yours is elsewhere

get_species_list() {
	case "$SET" in
	auto)
		# NOTE: depth is 4: ./<species>/<tier>/<inum>/polap-readassemble-1-miniasm
		find . -mindepth 4 -maxdepth 4 -type d \
			-path "./*/${TIER}/${INUM}/polap-readassemble-1-miniasm" -printf '%P\n' |
			cut -d/ -f1 | sort -u
		;;
	some)
		# Default "some" set (user-defined)
		cat <<'__SOME__'
Oryza_rufipogon
Trifolium_pratense
Dioscorea_japonica
Anthoceros_agrestis
Codonopsis_lanceolata
Canavalia_ensiformis
Arabidopsis_thaliana
Taraxacum_mongolicum
Cinchona_pubescens
Vitis_vinifera
Cucumis_sativus_var_hardwickii
Solanum_lycopersicum
Euonymus_alatus
Gossypium_herbaceum
Brassica_rapa
Phaeomegaceros_chiloensis
Juncus_effusus
Eucalyptus_pauciflora
Prunus_mandshurica
Juncus_inflexus
Juncus_roemerianus
Lolium_perenne
Dunaliella_tertiolecta
Notothylas_orbicularis
Anthoceros_angustus
Populus_x_sibirica
Spirodela_polyrhiza
Macadamia_jansenii
Vigna_radiata
Vaccinium_vitis_idaea
Carex_pseudochinensis
Leiosporoceros_dussii
Macadamia_tetraphylla
Musa_acuminata_subsp_malaccensis
Punica_granatum
Juncus_validus
Ophrys_lutea
Salix_dunnii
__SOME__
		;;
	file:*)
		local f="${SET#file:}"
		[[ -s "$f" ]] || {
			echo "[ERR] list not found: $f" >&2
			exit 2
		}
		sed -e 's/^[[:space:]]\+//' -e 's/[[:space:]]\+$//' -e '/^$/d' "$f"
		;;
	*)
		echo "$SET"
		;;
	esac
}

first_match() { compgen -G "$1" >/dev/null 2>&1 && ls $1 2>/dev/null | head -n1 || true; }

# ── WRITE-ONLY helper: echo one line to $FACTS (never fail) ───────────────
append_fact() {
	# species tier inum organelle kind key value
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$4" "$5" "$6" "$7" >>"$FACTS" || true
}

# ── SAFE per-species harvester (never exit nonzero) ──────────────────────────
harvest_one() {
	local sp="$1"
	{
		local base="${sp}/${TIER}/${INUM}/polap-readassemble-1-miniasm"
		# dataset stats (seqkit)
		local fq fqstat
		# Prefer precomputed stats if present (seqkit -T output)
		local fqstats="${sp}/${TIER}/${INUM}/summary-data/l.fq.stats"

		fq="$(first_match "${sp}/${TIER}/${INUM}/summary-data/l.fq.gz")"
		[[ -z "$fq" ]] && fq="$(first_match "${sp}/${TIER}/${INUM}/summary-data/l.fq")"
		if [[ -n "$fq" && -s "$fq" ]] && have seqkit; then
			fqstat="$(seqkit stats -T -a "$fq" 2>/dev/null | tail -n +2 || true)"
			if [[ -n "$fqstat" ]]; then
				IFS=$'\t' read -r _f _fmt _typ num sum _min avg _max n50 _q20 _q30 gc avgQ <<<"$fqstat"
				append_fact "$sp" "$TIER" "$INUM" "data" "attr" "total_bases" "${sum:-NA}"
				append_fact "$sp" "$TIER" "$INUM" "data" "attr" "read_count" "${num:-NA}"
				append_fact "$sp" "$TIER" "$INUM" "data" "attr" "mean_length" "${avg:-NA}"
				append_fact "$sp" "$TIER" "$INUM" "data" "attr" "N50" "${n50:-NA}"
				append_fact "$sp" "$TIER" "$INUM" "data" "attr" "avg_qual" "${avgQ:-NA}"
				append_fact "$sp" "$TIER" "$INUM" "data" "attr" "gc_content" "${gc:-NA}"
			fi

		elif [[ -s "$fqstats" ]]; then
			# Parse by header name, and emit TAB-separated fields
			local parsed
			parsed="$(
				awk -F'\t' '
          NR==1{
            for(i=1;i<=NF;i++) H[$i]=i
            next
          }
          NR==2{
            # Print exactly 6 fields in this order, TAB-separated:
            # num_seqs, sum_len, avg_len, N50, AvgQual, GC(%)
            printf "%s\t%s\t%s\t%s\t%s\t%s\n", \
              (H["num_seqs"]?$(H["num_seqs"]):"NA"), \
              (H["sum_len"]?$(H["sum_len"]):"NA"), \
              (H["avg_len"]?$(H["avg_len"]):"NA"), \
              (H["N50"]?$(H["N50"]):"NA"), \
              (H["AvgQual"]?$(H["AvgQual"]):"NA"), \
              (H["GC(%)"]?$(H["GC(%)"]):"NA")
          }
        ' "$fqstats"
			)"
			# Force TAB splitting here
			local num sum avg n50 avgQ gc
			IFS=$'\t' read -r num sum avg n50 avgQ gc <<<"$parsed"

			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "total_bases" "${sum:-NA}"
			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "read_count" "${num:-NA}"
			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "mean_length" "${avg:-NA}"
			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "N50" "${n50:-NA}"
			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "avg_qual" "${avgQ:-NA}"
			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "gc_content" "${gc:-NA}"
		fi

		# PT
		local pt_gfa
		pt_gfa="$(first_match "${base}/pt.1.gfa")"
		[[ -z "$pt_gfa" ]] && pt_gfa="$(first_match "${base}/pt1/assembly_graph.gfa")"
		if [[ -n "$pt_gfa" ]]; then
			append_fact "$sp" "$TIER" "$INUM" "pt" "file" "gfa" "$pt_gfa"
			append_fact "$sp" "$TIER" "$INUM" "pt" "file" "png" "${pt_gfa%.gfa}.png"
			local ctxt="${base}/pt1/ptdna/circular_path_count.txt"
			[[ -s "$ctxt" ]] && append_fact "$sp" "$TIER" "$INUM" "pt" "attr" "circular_count" "$(cat "$ctxt" 2>/dev/null || echo NA)"
			local cfa
			cfa="$(first_match "${base}/pt1/ptdna/circular_path_*_concatenated.fa")"
			[[ -n "$cfa" ]] && append_fact "$sp" "$TIER" "$INUM" "pt" "file" "circular_fa" "$cfa"
			local gc
			gc="$(first_match "${base}/pt1/50-annotation/pt.gene.count")"
			[[ -n "$gc" ]] && append_fact "$sp" "$TIER" "$INUM" "pt" "attr" "gene_count" "$(cat "$gc" 2>/dev/null || echo NA)" &&
				append_fact "$sp" "$TIER" "$INUM" "pt" "file" "gene_count_file" "$gc"
			local sc
			sc="$(first_match "${base}/annotate-read-mtseed/pt/pt-contig-annotation-depth-table.txt.scatter.pdf")"
			[[ -z "$sc" ]] && sc="$(first_match "${base}/annotate-read-pt/pt/pt-contig-annotation-depth-table.txt.scatter.pdf")"
			[[ -n "$sc" ]] && append_fact "$sp" "$TIER" "$INUM" "pt" "file" "qc_scatter" "$sc"
		fi

		# ── PT ────────────────────────────────────────────────────────────────────────
		# try in order: pt.1.gfa, pt-pt.1.gfa, pt1/assembly_graph.gfa
		pt_gfa="$(first_match "${base}/pt.1.gfa")"
		[[ -z "$pt_gfa" ]] && pt_gfa="$(first_match "${base}/pt-pt.1.gfa")"
		[[ -z "$pt_gfa" ]] && pt_gfa="$(first_match "${base}/pt1/assembly_graph.gfa")"

		if [[ -n "$pt_gfa" ]]; then
			pt_png="${pt_gfa%.gfa}.png"
			# render PNG if missing (Bandage via polap.sh)
			if [[ ! -s "$pt_png" && -s "$pt_gfa" ]] && command -v polap.sh >/dev/null 2>&1; then
				polap.sh bandage png "$pt_gfa" "$pt_png" || true
			fi
			append_fact "$sp" "$TIER" "$INUM" "pt" "file" "gfa" "$pt_gfa"
			append_fact "$sp" "$TIER" "$INUM" "pt" "file" "png" "$pt_png"
		fi

		# PT-PT
		if [[ $INCLUDE_PTPT -eq 1 ]]; then
			local ptpt
			ptpt="$(first_match "${base}/pt-pt.1.gfa")"
			[[ -n "$ptpt" ]] && append_fact "$sp" "$TIER" "$INUM" "ptpt" "file" "gfa" "$ptpt" &&
				append_fact "$sp" "$TIER" "$INUM" "ptpt" "file" "png" "${ptpt%.gfa}.png"
		fi

		# MT
		local mt_gfa
		mt_gfa="$(first_match "${base}/mt.1.gfa")"
		[[ -z "$mt_gfa" ]] && mt_gfa="$(first_match "${base}/mt1/assembly_graph.gfa")"
		if [[ -n "$mt_gfa" ]]; then
			append_fact "$sp" "$TIER" "$INUM" "mt" "file" "gfa" "$mt_gfa"
			append_fact "$sp" "$TIER" "$INUM" "mt" "file" "png" "${mt_gfa%.gfa}.png"
			local mgc
			mgc="$(first_match "${base}/mt1/50-annotation/mt.gene.count")"
			[[ -n "$mgc" ]] && append_fact "$sp" "$TIER" "$INUM" "mt" "attr" "gene_count" "$(cat "$mgc" 2>/dev/null || echo NA)" &&
				append_fact "$sp" "$TIER" "$INUM" "mt" "file" "gene_count_file" "$mgc"
			local msc
			msc="$(first_match "${base}/annotate-read-mtseed/mt/contig-annotation-depth-table.txt.scatter.pdf")"
			[[ -n "$msc" ]] && append_fact "$sp" "$TIER" "$INUM" "mt" "file" "qc_scatter" "$msc"
		fi

	} || true # NEVER kill the whole job
}

# ── parallel-safe fanout: per-species files, then merge ──────────────────────
FACTS="$(mktemp -t polap-facts.XXXXXX.tsv)"
echo -e "species\ttier\tinum\torganelle\tkind\tkey\tvalue" >"$FACTS"
FACTS_DIR="$(mktemp -d -t polap-factsdir.XXXXXX)"
trap 'rm -f "$FACTS"; rm -rf "$FACTS_DIR"' EXIT

harvest_wrapper() {
	local sp="$1" out="$FACTS_DIR/${sp}.tsv"
	FACTS="$out" harvest_one "$sp"
}
export -f harvest_one harvest_wrapper first_match append_fact
export TIER INUM INCLUDE_PTPT _POLAPLIB_DIR

# Try parallel; fallback to serial
# if have parallel; then
# 	get_species_list | parallel --will-cite -j"$(nproc)" --halt soon,fail=1 bash -c 'harvest_wrapper "$@"' _ {} || true
# else
while read -r sp; do harvest_wrapper "$sp"; done < <(get_species_list)
# fi

# Merge per-species
if ls "$FACTS_DIR"/*.tsv >/dev/null 2>&1; then LC_ALL=C sort -V "$FACTS_DIR"/*.tsv >>"$FACTS"; fi

# Assemble JSON
python3 "${_POLAPLIB_DIR}/scripts/manifest_assemble.py" \
	--facts "$FACTS" \
	--set "$SET" --tier "$TIER" --inum "$INUM" \
	--out "$OUT" \
	$([[ $PRETTY -eq 1 ]] && echo --pretty)

echo "[INFO] Wrote: $OUT"
