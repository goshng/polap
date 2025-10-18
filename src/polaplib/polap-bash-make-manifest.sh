#!/usr/bin/env bash
# polaplib/polap-bash-make-manifest.sh
#
# Version : v0.4.2  (2025-10-15)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Purpose:
#   Harvest per-species facts and assemble a structured JSON manifest via
#   scripts/manifest_assemble.py.
#
# This version (v0.4.2):
#   • Writes long-read block `data` and short-read blocks `short1data`/`short2data`
#     even when NO seqkit stats exist — by emitting placeholder rows.
#     (at minimum: sra_id="" and total_bases="")
#   • Keeps robust parsing of seqkit -Ta tables and SRA detection.
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
QUIET="${QUIET:-1}"

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --set some|auto|file:/path|Species_Name
                   --tier v6 --inum 0
                   --out md/manifest.json
                   [--include-ptpt] [--pretty] [--debug]
EOF
}

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
		find . -mindepth 4 -maxdepth 4 -type d \
			-path "./*/${TIER}/${INUM}/polap-readassemble" -printf '%P\n' |
			cut -d/ -f1 | LC_ALL=C sort -u
		;;
	some)
		cat <<'__SOME__'
Anthoceros_agrestis
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

# Parse a single “seqkit stats -T/-Ta” table → prints: num  sum  avg  n50  avgQ  gc
parse_seqkit_ta() {
	local f="$1"
	awk '
	  function norm(s){gsub(/\r$/,"",s); gsub(/^[[:space:]]+|[[:space:]]+$/,"",s); return tolower(s)}
	  BEGIN{FS="\t"}
	  function refit(){
	    if (NF==1){ FS="[ \t]+"; n=split($0,a,FS); for(i=1;i<=n;i++) $i=a[i]; NF=n }
	  }
	  NR==1{
	    refit()
	    for(i=1;i<=NF;i++){ key=norm($i); H[key]=i; if(key ~ /^gc\(%\)$/) GC=i }
	    next
	  }
	  NR==2{
	    refit()
	 	n   = (H["num_seqs"]?$(H["num_seqs"]):"")
	    s   = (H["sum_len"] ?$(H["sum_len"]) :"")
	    a   = (H["avg_len"] ?$(H["avg_len"]) :"")
	    n50 = (H["n50"]     ?$(H["n50"])     :"")
	    aq  = (H["avgqual"] ?$(H["avgqual"]):(H["avg_qual"]?$(H["avg_qual"]):""))
	    if (GC) gc=$(GC); else if(H["gc(%)"]) gc=$(H["gc(%)"]); else if(H["gc"]) gc=$(H["gc"]); else gc=""
	    printf "%s\t%s\t%s\t%s\t%s\t%s\n", n,s,a,n50,aq,gc; exit
	  }
	' "$f"
}

# Fallback: compute stats one-liner (quiet)
compute_seqkit_Ta() {
	local in="$1"
	seqkit stats -Ta "$in" 2>/dev/null | tail -n +2 || true
}

_emit_six_attrs() { # species tier inum tag num sum avg n50 avgQ gc
	local sp="$1" tier="$2" inum="$3" tag="$4"
	local num="$5" sum="$6" avg="$7" n50="$8" avgQ="$9" gc="${10}"
	append_fact "$sp" "$tier" "$inum" "$tag" "attr" "total_bases" "${sum:-}"
	append_fact "$sp" "$tier" "$inum" "$tag" "attr" "read_count" "${num:-}"
	append_fact "$sp" "$tier" "$inum" "$tag" "attr" "mean_length" "${avg:-}"
	append_fact "$sp" "$tier" "$inum" "$tag" "attr" "N50" "${n50:-}"
	append_fact "$sp" "$tier" "$inum" "$tag" "attr" "avg_qual" "${avgQ:-}"
	append_fact "$sp" "$tier" "$inum" "$tag" "attr" "gc_content" "${gc:-}"
}

# Ensure a block exists even when no data: append empty sra_id and total_bases
_ensure_empty_block() {
	local sp="$1" tier="$2" inum="$3" tag="$4"
	append_fact "$sp" "$tier" "$inum" "$tag" "attr" "sra_id" ""
	append_fact "$sp" "$tier" "$inum" "$tag" "attr" "total_bases" ""
}

# ------------------------------------------------------------------------------
# Harvest one species (serial; never abort the run)
# ------------------------------------------------------------------------------
harvest_one() {
	local sp="$1"
	{
		local base="${sp}/${TIER}/${INUM}/polap-readassemble"
		local sdir="${sp}/${TIER}/${INUM}/summary-data"

		# Track if blocks got any rows
		local have_data=0 have_s1=0 have_s2=0

		# ---------- LONG: metrics from l.fq.seqkit.stats.ta.txt (or fallback) ---
		local l_ta="${sdir}/l.fq.seqkit.stats.ta.txt"
		if [[ -s "$l_ta" ]]; then
			local row
			row="$(parse_seqkit_ta "$l_ta")"
			if [[ -n "$row" ]]; then
				IFS=$'\t' read -r num sum avg n50 avgQ gc <<<"$row"
				_emit_six_attrs "$sp" "$TIER" "$INUM" "data" "$num" "$sum" "$avg" "$n50" "$avgQ" "$gc"
				have_data=1
			fi
		else
			local fq="$(first_match "${sdir}/l.fq.gz")"
			[[ -z "$fq" ]] && fq="$(first_match "${sdir}/l.fq")"
			if [[ -n "$fq" && -s "$fq" ]] && have seqkit; then
				local row
				row="$(compute_seqkit_Ta "$fq")"
				if [[ -n "$row" ]]; then
					IFS=$'\t' read -r _f _fmt _typ num sum _min avg _max _Q1 _Q2 _Q3 _sumgap n50 _n50num _q20 _q30 avgQ gc _sumn <<<"$row"
					_emit_six_attrs "$sp" "$TIER" "$INUM" "data" "$num" "$sum" "$avg" "$n50" "$avgQ" "$gc"
					have_data=1
				fi
			fi
		fi

		# data_file_bytes from <species>/tmp/l.fq (uncompressed, if exists)
		local lf="$sp/tmp/l.fq"
		if [[ -s "$lf" ]]; then
			local bytes=""
			bytes="$(wc -c <"$lf" | tr -d '[:space:]' || true)"
			[[ -n "$bytes" ]] && append_fact "$sp" "$TIER" "$INUM" "data" "attr" "data_file_bytes" "$bytes" && have_data=1
		fi

		# LONG SRA: prefer tmp/l.sra.txt, fallback summary-data
		local l_sra=""
		[[ -s "$sp/tmp/l.sra.txt" ]] && l_sra="$(head -n1 "$sp/tmp/l.sra.txt" | tr -d '\r\n')"
		[[ -z "$l_sra" && -s "$sdir/l.sra.txt" ]] && l_sra="$(head -n1 "$sdir/l.sra.txt" | tr -d '\r\n')"
		if [[ -n "$l_sra" ]]; then
			append_fact "$sp" "$TIER" "$INUM" "data" "attr" "sra_id" "$l_sra"
			have_data=1
		fi

		# ---------- SHORT SRA (pair): prefer tmp/s.sra.txt, fallback summary-data
		local s_sra=""
		[[ -s "$sp/tmp/s.sra.txt" ]] && s_sra="$(head -n1 "$sp/tmp/s.sra.txt" | tr -d '\r\n')"
		[[ -z "$s_sra" && -s "$sdir/s.sra.txt" ]] && s_sra="$(head -n1 "$sdir/s.sra.txt" | tr -d '\r\n')"
		if [[ -n "$s_sra" ]]; then
			append_fact "$sp" "$TIER" "$INUM" "short1data" "attr" "sra_id" "$s_sra"
			have_s1=1
			append_fact "$sp" "$TIER" "$INUM" "short2data" "attr" "sra_id" "$s_sra"
			have_s2=1
		fi

		# ---------- SHORT metrics for mates -------------------------------------
		for mate in 1 2; do
			local tag="short${mate}data"
			local s_ta="${sdir}/s_${mate}.fq.seqkit.stats.ta.txt"
			if [[ -s "$s_ta" ]]; then
				local row
				row="$(parse_seqkit_ta "$s_ta")"
				if [[ -n "$row" ]]; then
					IFS=$'\t' read -r num sum avg n50 avgQ gc <<<"$row"
					_emit_six_attrs "$sp" "$TIER" "$INUM" "$tag" "$num" "$sum" "$avg" "$n50" "$avgQ" "$gc"
					[[ "$mate" == "1" ]] && have_s1=1 || have_s2=1
				fi
			else
				local sfq="$(first_match "${sdir}/s_${mate}.fq.gz")"
				[[ -z "$sfq" ]] && sfq="$(first_match "${sdir}/s_${mate}.fq")"
				[[ -z "$sfq" ]] && sfq="$(first_match "${sp}/tmp/s_${mate}.fq.gz")"
				[[ -z "$sfq" ]] && sfq="$(first_match "${sp}/tmp/s_${mate}.fq")"
				if [[ -n "$sfq" && -s "$sfq" ]] && have seqkit; then
					local row
					row="$(compute_seqkit_Ta "$sfq")"
					if [[ -n "$row" ]]; then
						IFS=$'\t' read -r _f _fmt _typ num sum _min avg _max _Q1 _Q2 _Q3 _sumgap n50 _n50num _q20 _q30 avgQ gc _sumn <<<"$row"
						_emit_six_attrs "$sp" "$TIER" "$INUM" "$tag" "$num" "$sum" "$avg" "$n50" "$avgQ" "$gc"
						[[ "$mate" == "1" ]] && have_s1=1 || have_s2=1
					fi
				fi
			fi
			# If we already had s_sra, ensure it’s present on each short block
			if [[ -n "$s_sra" ]]; then
				append_fact "$sp" "$TIER" "$INUM" "$tag" "attr" "sra_id" "$s_sra"
				[[ "$mate" == "1" ]] && have_s1=1 || have_s2=1
			fi
		done

		# ---------- Ensure empty entries for missing blocks ----------------------
		((have_data == 0)) && _ensure_empty_block "$sp" "$TIER" "$INUM" "data"
		((have_s1 == 0)) && _ensure_empty_block "$sp" "$TIER" "$INUM" "short1data"
		((have_s2 == 0)) && _ensure_empty_block "$sp" "$TIER" "$INUM" "short2data"

		# ------------------------------ PT files ---------------------------------
		local base_asm="${sp}/${TIER}/${INUM}/polap-readassemble"
		local pt_gfa
		pt_gfa="$(first_match "${base_asm}/pt.1.gfa")"
		[[ -z "$pt_gfa" ]] && pt_gfa="$(first_match "${base_asm}/pt-pt.1.gfa")"
		[[ -z "$pt_gfa" ]] && pt_gfa="$(first_match "${base_asm}/pt1/assembly_graph.gfa")"
		if [[ -n "$pt_gfa" ]]; then
			append_fact "$sp" "$TIER" "$INUM" "pt" "file" "gfa" "$pt_gfa"
			append_fact "$sp" "$TIER" "$INUM" "pt" "file" "png" "${pt_gfa%.gfa}.png"
			local ctxt="${base_asm}/pt1/ptdna/circular_path_count.txt"
			[[ -s "$ctxt" ]] && append_fact "$sp" "$TIER" "$INUM" "pt" "attr" "circular_count" "$(cat "$ctxt" 2>/dev/null || echo "")"
			local cfa
			cfa="$(first_match "${base_asm}/pt1/ptdna/circular_path_*_concatenated.fa")"
			[[ -n "$cfa" ]] && append_fact "$sp" "$TIER" "$INUM" "pt" "file" "circular_fa" "$cfa"
			local gc
			gc="$(first_match "${base_asm}/pt1/50-annotation/pt.gene.count")"
			[[ -n "$gc" ]] && append_fact "$sp" "$TIER" "$INUM" "pt" "attr" "gene_count" "$(cat "$gc" 2>/dev/null || echo "")" &&
				append_fact "$sp" "$TIER" "$INUM" "pt" "file" "gene_count_file" "$gc"
			local sc
			sc="$(first_match "${base_asm}/annotate-read-mtseed/pt/pt-contig-annotation-depth-table.txt.scatter.pdf")"
			[[ -z "$sc" ]] && sc="$(first_match "${base_asm}/annotate-read-pt/pt/pt-contig-annotation-depth-table.txt.scatter.pdf")"
			[[ -n "$sc" ]] && append_fact "$sp" "$TIER" "$INUM" "pt" "file" "qc_scatter" "$sc"
		fi

		# PT ref accession/length (typo tolerant)
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

		# ------------------------------ MT files ---------------------------------
		local mt_gfa
		mt_gfa="$(first_match "${base_asm}/mt.1.gfa")"
		[[ -z "$mt_gfa" ]] && mt_gfa="$(first_match "${base_asm}/mt1/assembly_graph.gfa")"
		if [[ -n "$mt_gfa" ]]; then
			append_fact "$sp" "$TIER" "$INUM" "mt" "file" "gfa" "$mt_gfa"
			append_fact "$sp" "$TIER" "$INUM" "mt" "file" "png" "${mt_gfa%.gfa}.png"
			local mgc
			mgc="$(first_match "${base_asm}/mt1/50-annotation/mt.gene.count")"
			[[ -n "$mgc" ]] && append_fact "$sp" "$TIER" "$INUM" "mt" "attr" "gene_count" "$(cat "$mgc" 2>/dev/null || echo "")" &&
				append_fact "$sp" "$TIER" "$INUM" "mt" "file" "gene_count_file" "$mgc"
			local msc
			msc="$(first_match "${base_asm}/annotate-read-mtseed/mt/contig-annotation-depth-table.txt.scatter.pdf")"
			[[ -n "$msc" ]] && append_fact "$sp" "$TIER" "$INUM" "mt" "file" "qc_scatter" "$msc"
		fi

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

while read -r sp; do
	[[ -z "${sp// /}" ]] && continue
	[[ $DEBUG -eq 1 ]] && echo "[DBG] harvesting: $sp" >&2
	harvest_one "$sp" >/dev/null
done < <(get_species_list)

[[ $DEBUG -eq 1 ]] && {
	echo "[DBG] FACTS at: $FACTS"
	grep -E $'\t(data|short1data|short2data)\tattr\t(sra_id|total_bases)' "$FACTS" || true
}

# assemble
CODES_FILE="${_POLAPLIB_DIR}/species-codes.txt"
python3 "${_POLAPLIB_DIR}/scripts/manifest_assemble.py" \
	--facts "$FACTS" \
	--set "$SET" --tier "$TIER" --inum "$INUM" \
	--out "$OUT" $([[ $PRETTY -eq 1 ]] && echo --pretty) \
	--codes "$CODES_FILE"

[[ "$QUIET" -eq 1 ]] || echo "[INFO] Wrote: $OUT"
