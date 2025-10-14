#!/usr/bin/env bash
# polaplib/polap-bash-make-manifest.sh
#
# Version : v0.3.0  (2025-10-13)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Purpose:
#   Harvest per-species facts (dataset stats, PT/MT files, NCBI PT ref)
#   and assemble a structured JSON manifest via scripts/manifest_assemble.py.
#
# Notes:
#   - Serial harvesting (no GNU parallel)
#   - Header-aware seqkit -T parsing (num_seqs,sum_len,avg_len,N50,AvgQual,GC(%))
#   - Records data.sra_id, data.data_file_bytes (species/tmp/l.fq)
#   - PT/MT GFA, PNG, ancillary files (qc/gene/circular)
#   - Adds PT ref: ncbi_accession, ncbi_len_bp from ncbi-ptdna/00-bioproject
#   - Tolerates “mtdna” typo in PT ref stats file names
#
set -euo pipefail
IFS=$'\n\t'

# ------------------------------------------------------------------------------
# Auto base & env
# ------------------------------------------------------------------------------
_POLAPLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export _POLAPLIB_DIR

# ------------------------------------------------------------------------------
# CLI
# ------------------------------------------------------------------------------
SET="some" # auto|some|file:/path|<Species_Name>
TIER="v6"
INUM="0"
OUT="md/manifest.json"
INCLUDE_PTPT=0
PRETTY=0
DEBUG=0

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --set some|auto|file:/path|Species_Name
                   --tier v6 --inum 0
                   --out md/manifest.json
                   [--include-ptpt] [--pretty] [--debug]
EOF
}

# at top defaults
QUIET="${QUIET:-1}"

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
	--debug)
		DEBUG=1
		shift
		;;
	-h | --help)
		usage
		exit 0
		;;
	*)
		echo "[ERR] Unknown arg: $1" >&2
		usage
		exit 2
		;;
	esac
done
mkdir -p "$(dirname "$OUT")"

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
have() { command -v "$1" >/dev/null 2>&1; }
first_match() { compgen -G "$1" >/dev/null 2>&1 && ls $1 2>/dev/null | head -n1 || true; }

get_species_list() {
	case "$SET" in
	auto)
		# depth=4  ./<species>/<tier>/<inum>/polap-readassemble-1-miniasm
		find . -mindepth 4 -maxdepth 4 -type d \
			-path "./*/${TIER}/${INUM}/polap-readassemble-1-miniasm" -printf '%P\n' |
			cut -d/ -f1 | LC_ALL=C sort -u
		;;
	some)
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
	*) echo "$SET" ;;
	esac
}

append_fact() {
	# species tier inum organelle kind key value
	printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$4" "$5" "$6" "$7" >>"$FACTS" || true
}

# ------------------------------------------------------------------------------
# Harvest one species (serial; never abort the run)
# ------------------------------------------------------------------------------
harvest_one() {
	local sp="$1"
	{
		local base="${sp}/${TIER}/${INUM}/polap-readassemble-1-miniasm"

		# ---------- dataset stats (seqkit or precomputed l.fq.stats) -------------
		local fq fqstats fqstat
		fqstats="${sp}/${TIER}/${INUM}/summary-data/l.fq.stats"

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
			local parsed num sum avg n50 avgQ gc
			parsed="$(
				awk -F'\t' '
          NR==1{for(i=1;i<=NF;i++)H[$i]=i; next}
          NR==2{printf "%s\t%s\t%s\t%s\t%s\t%s\n", \
            (H["num_seqs"]?$(H["num_seqs"]):"NA"), \
            (H["sum_len"]? $(H["sum_len"]) :"NA"), \
            (H["avg_len"]? $(H["avg_len"]) :"NA"), \
            (H["N50"]?    $(H["N50"])    :"NA"), \
            (H["AvgQual"]?$(H["AvgQual"]):"NA"), \
            (H["GC(%)"]?  $(H["GC(%)"])  :"NA") }' "$fqstats"
			)"
			IFS=$'\t' read -r num sum avg n50 avgQ gc <<<"$parsed"
			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "total_bases" "${sum:-NA}"
			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "read_count" "${num:-NA}"
			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "mean_length" "${avg:-NA}"
			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "N50" "${n50:-NA}"
			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "avg_qual" "${avgQ:-NA}"
			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "gc_content" "${gc:-NA}"
		fi

		# SRA id
		local sra_txt="${sp}/${TIER}/${INUM}/summary-data/l.sra.txt"
		if [[ -s "$sra_txt" ]]; then
			local sra_id
			sra_id="$(head -n1 "$sra_txt" | tr -d '\r\n')"
			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "sra_id" "${sra_id:-NA}"
		fi

		# data bytes (species/tmp/l.fq)
		local lf="$sp/tmp/l.fq"
		if [[ -s "$lf" ]]; then
			local bytes="NA"
			if stat --printf='%s' "$lf" >/dev/null 2>&1; then
				bytes="$(stat --printf='%s' "$lf")"
			else
				bytes="$(wc -c <"$lf" | tr -d '[:space:]')"
			fi
			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "data_file_bytes" "${bytes}"
		fi

		# ------------------------------ PT files ---------------------------------
		local pt_gfa
		pt_gfa="$(first_match "${base}/pt.1.gfa")"
		[[ -z "$pt_gfa" ]] && pt_gfa="$(first_match "${base}/pt-pt.1.gfa")"
		[[ -z "$pt_gfa" ]] && pt_gfa="$(first_match "${base}/pt1/assembly_graph.gfa")"
		if [[ -n "$pt_gfa" ]]; then
			local pt_png="${pt_gfa%.gfa}.png"
			if [[ ! -s "$pt_png" && -s "$pt_gfa" ]] && have polap.sh; then
				polap.sh bandage png "$pt_gfa" "$pt_png" || true
			fi
			append_fact "$sp" "$TIER" "$INUM" "pt" "file" "gfa" "$pt_gfa"
			append_fact "$sp" "$TIER" "$INUM" "pt" "file" "png" "$pt_png"

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

		# PT ref accession/length from ncbi-ptdna/00-bioproject (typo tolerant)
		local biop="${sp}/${TIER}/${INUM}/ncbi-ptdna/00-bioproject"
		if [[ -d "$biop" ]]; then
			local stats
			stats="$(first_match "${biop}/1-ptdna.fasta.stats")"
			[[ -z "$stats" ]] && stats="$(first_match "${biop}/1-mtdna.fasta.stats")"
			[[ -z "$stats" ]] && stats="$(first_match "${biop}/ptdna.fasta.stats")"
			[[ -z "$stats" ]] && stats="$(first_match "${biop}"/*fasta.stats)"
			if [[ -n "$stats" ]]; then
				local acc len
				acc="$(awk 'NR==2{print $1; exit}' "$stats")"
				len="$(awk 'NR==2{print $NF; exit}' "$stats")"
				[[ -n "$acc" ]] && append_fact "$sp" "$TIER" "$INUM" "pt" "attr" "ncbi_accession" "$acc"
				[[ -n "$len" ]] && append_fact "$sp" "$TIER" "$INUM" "pt" "attr" "ncbi_len_bp" "$len"
			else
				local accf
				accf="$(first_match "${biop}/2-ptdna.accession")"
				[[ -z "$accf" ]] && accf="$(first_match "${biop}/2-mtdna.accession")"
				[[ -z "$accf" ]] && accf="$(first_match "${biop}"/*.accession)"
				if [[ -n "$accf" ]]; then
					local acc
					acc="$(head -n1 "$accf" | tr -d '\r\n')"
					[[ -n "$acc" ]] && append_fact "$sp" "$TIER" "$INUM" "pt" "attr" "ncbi_accession" "$acc"
				fi
			fi
		fi

		if [[ $INCLUDE_PTPT -eq 1 ]]; then
			local ptpt
			ptpt="$(first_match "${base}/pt-pt.1.gfa")"
			[[ -n "$ptpt" ]] && append_fact "$sp" "$TIER" "$INUM" "ptpt" "file" "gfa" "$ptpt" &&
				append_fact "$sp" "$TIER" "$INUM" "ptpt" "file" "png" "${ptpt%.gfa}.png"
		fi

		# ------------------------------ MT files ---------------------------------
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

		# NCBI MT accession/length from ncbi-mtdna/00-bioproject
		local biomt="${sp}/${TIER}/${INUM}/ncbi-mtdna/00-bioproject"
		if [[ -d "$biomt" ]]; then
			local mstats
			mstats="$(first_match "${biomt}/1-mtdna.fasta.stats")"
			[[ -z "$mstats" ]] && mstats="$(first_match "${biomt}/mtdna.fasta.stats")"
			[[ -z "$mstats" ]] && mstats="$(first_match "${biomt}"/*fasta.stats)"
			if [[ -n "$mstats" ]]; then
				local macc mlen
				macc="$(awk 'NR==2{print $1; exit}' "$mstats")"
				mlen="$(awk 'NR==2{print $NF; exit}' "$mstats")"
				[[ -n "$macc" ]] && append_fact "$sp" "$TIER" "$INUM" "mt" "attr" "ncbi_accession" "$macc"
				[[ -n "$mlen" ]] && append_fact "$sp" "$TIER" "$INUM" "mt" "attr" "ncbi_len_bp" "$mlen"
			else
				local maccf
				maccf="$(first_match "${biomt}/2-mtdna.accession")"
				[[ -z "$maccf" ]] && maccf="$(first_match "${biomt}"/*.accession)"
				if [[ -n "$maccf" ]]; then
					local macc
					macc="$(head -n1 "$maccf" | tr -d '\r\n')"
					[[ -n "$macc" ]] && append_fact "$sp" "$TIER" "$INUM" "mt" "attr" "ncbi_accession" "$macc"
				fi
			fi
		fi

	} || true
}

# ------------------------------------------------------------------------------
# Collect facts (serial), then assemble JSON
# ------------------------------------------------------------------------------
FACTS="$(mktemp -t polap-facts.XXXXXX.tsv)"
echo -e "species\ttier\tinum\torganelle\tkind\tkey\tvalue" >"$FACTS"

# serial harvest
while read -r sp; do
	[[ -z "${sp// /}" ]] && continue
	[[ $DEBUG -eq 1 ]] && echo "[DBG] harvesting: $sp" >&2
	harvest_one "$sp" >/dev/null
done < <(get_species_list)

[[ $DEBUG -eq 1 ]] && {
	echo "[DBG] FACTS at: $FACTS"
	head -n 25 "$FACTS" >&2
}

# assemble
CODES_FILE="${_POLAPLIB_DIR}/species-codes.txt"

# assemble JSON
python3 "${_POLAPLIB_DIR}/scripts/manifest_assemble.py" \
	--facts "$FACTS" \
	--set "$SET" --tier "$TIER" --inum "$INUM" \
	--out "$OUT" $([[ $PRETTY -eq 1 ]] && echo --pretty) \
	--codes "$CODES_FILE"
# $([[ -s "$CODES_FILE" ]] && echo --codes "$CODES_FILE")

# only print if not in quiet mode
[[ "$QUIET" -eq 1 ]] || echo "[INFO] Wrote: $OUT"
