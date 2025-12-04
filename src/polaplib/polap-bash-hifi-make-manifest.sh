#!/usr/bin/env bash
# a copy of polaplib/polap-bash-make-manifest.sh
#
# polaplib/polap-bash-hifi-make-manifest.sh
#
# Version : v0.5.0  (2025-10-15)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Purpose:
#   Harvest per-species facts and assemble a structured JSON manifest via
#   scripts/manifest_assemble.py.
#
# This version (v0.5.0):
#   • Writes long-read block `data` and short-read blocks `short1data`/`short2data`
#     even when NO seqkit stats exist — by emitting placeholder rows.
#     (at minimum: sra_id="" and total_bases="")
#   • Adds HiFi assembly entries for OATK (all summary-oatk-*.txt variants),
#     TIPPo, HiMT, and PMAT2, including:
#       - mtDNA/ptDNA GFA, FASTA, PNG paths
#       - resource usage (mem/time/disk) from summary-*.txt
#       - HiMT 4C metrics (geneset completeness, contiguity, fragmentation index,
#         largest contig proportion, contig-length CV) parsed from assess/<tool>/report/prefix.* files.
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
# This codes should be from the manifest.
SPECIES_CODES_TXT="polap-species-codes-read.txt"
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
                   --species-codes FILE
                   --out md/manifest.json
                   [--include-ptpt] [--pretty] [--debug]
EOF
}

while (($#)); do
	case "$1" in
	--species-codes)
		SPECIES_CODES_TXT="${2:?}"
		shift 2
		;;
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

# species list (HiFi/root)
get_species_list() {
	case "$SET" in
	auto)
		find . -mindepth 4 -maxdepth 4 -type d \
			-path "./*/${TIER}/${INUM}/oatk" -printf '%P\n' |
			cut -d/ -f1 | LC_ALL=C sort -u
		;;
	autoread)
		find . -mindepth 4 -maxdepth 4 -type d \
			-path "./*/${TIER}/${INUM}/polap-readassemble" -printf '%P\n' |
			cut -d/ -f1 | LC_ALL=C sort -u
		;;
	some)
		cat <<'__SOME__'
Biancaea_sappan
__SOME__
		;;
	someread)
		cat <<'__SOME__'
Anthoceros_agrestis
Anthoceros_angustus
Arabidopsis_thaliana
Codonopsis_lanceolata
Eucalyptus_pauciflora
Juncus_roemerianus
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

# Parse a single “seqkit stats -T/-Ta” table -> prints: num  sum  avg  n50  avgQ  gc
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

# Parse summary-*.txt: set globals sum_mem_kb, sum_mem_gb, sum_time_hms, sum_time_hours, sum_disk_gb
_extract_summary_metrics() {
	local f="$1"
	sum_mem_kb=""
	sum_mem_gb=""
	sum_time_hms=""
	sum_time_hours=""
	sum_disk_gb=""

	[[ ! -s "$f" ]] && return 0

	local line

	# Net increase memory line, e.g.:
	# Net increase:    9316600 KB (8.89 GB)
	line="$(grep -m1 -E 'Net increase' "$f" || true)"
	if [[ -n "$line" ]]; then
		sum_mem_kb="$(awk '{for(i=1;i<=NF;i++) if($i=="KB"||$i=="KB)"){print $(i-1);exit}}' <<<"$line" | tr -d '[:space:]')"
		sum_mem_gb="$(sed -n 's/.*(\(.*GB\)).*/\1/p' <<<"$line" | tr -d '[:space:]')"
	fi

	# Elapsed time line, e.g.:
	# Elapsed time:    00:23:42 (0.40 h)
	line="$(grep -m1 -E 'Elapsed time' "$f" || true)"
	if [[ -n "$line" ]]; then
		sum_time_hms="$(awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]{2}:[0-9]{2}:[0-9]{2}$/){print $i;exit}}' <<<"$line")"
		sum_time_hours="$(sed -n 's/.*(\(.*h\)).*/\1/p' <<<"$line" | tr -d '[:space:]')"
		# Fallback if no parentheses
		if [[ -z "$sum_time_hours" && -n "$sum_time_hms" ]]; then
			IFS=: read -r hh mm ss <<<"$sum_time_hms"
			local total_sec=$((10#$hh * 3600 + 10#$mm * 60 + 10#$ss))
			sum_time_hours="$(awk -v sec="$total_sec" 'BEGIN{printf"%.2fh", sec/3600}')"
		fi
	fi

	# Disk used line, e.g.:
	#   Disk used:       24 GB
	#   Disk used:       24GB
	#   Disk used:       24.5 GB
	line="$(grep -m1 -E '^[[:space:]]*Disk[[:space:]]+used' "$f" || true)"
	sum_disk_gb=""
	if [[ -n "$line" ]]; then
		# Extract the first number that is followed by optional space + GB
		sum_disk_gb="$(
			printf '%s\n' "$line" |
				sed -E 's/.*Disk[[:space:]]+used:[^0-9]*([0-9]+(\.[0-9]+)?)[[:space:]]*GB.*/\1/'
		)"
		# If that didn't work (e.g. 24GB with no space), try a more generic fallback:
		if [[ -z "$sum_disk_gb" ]]; then
			sum_disk_gb="$(
				printf '%s\n' "$line" |
					sed -E 's/.*Disk[[:space:]]+used:[^0-9]*([0-9]+(\.[0-9]+)?)GB.*/\1/'
			)"
		fi
	fi

}

# Parse HiMT metrics in assess/<tool>/report/prefix.*; set globals:
#   hmt_geneset, hmt_contigs, hmt_total_len, hmt_n50, hmt_frag_idx,
#   hmt_max_prop, hmt_cv
_extract_himt_metrics() {
	local sp="$1" tier="$2" inum="$3" tool="$4"
	hmt_geneset=""
	hmt_contigs=""
	hmt_total_len=""
	hmt_n50=""
	hmt_frag_idx=""
	hmt_max_prop=""
	hmt_cv=""

	local rdir="${sp}/${tier}/${inum}/assess/${tool}/report"
	local gis="${rdir}/prefix.gene_integrity_summary.tsv"
	local mb="${rdir}/prefix.mito_basic_info.tsv"
	local ct="${rdir}/prefix.contig_table.tsv"

	# geneset completeness prop
	if [[ -s "$gis" ]]; then
		hmt_geneset="$(awk -F'\t' '$1=="geneset_completeness_prop"{print $2}' "$gis" | head -n1)"
	fi

	# mito basic info: total length, contig number, N50, min/max length
	local maxlen="" minlen=""
	if [[ -s "$mb" ]]; then
		hmt_total_len="$(awk -F'\t' '$1=="Total length"{print $2}' "$mb" | head -n1)"
		hmt_contigs="$(awk -F'\t' '$1=="Total contig number"{print $2}' "$mb" | head -n1)"
		hmt_n50="$(awk -F'\t' '$1=="N50"{print $2}' "$mb" | head -n1)"
		maxlen="$(awk -F'\t' '$1=="Maximum length"{print $2}' "$mb" | head -n1)"
		minlen="$(awk -F'\t' '$1=="Minimum length"{print $2}' "$mb" | head -n1)"
	fi

	# fragmentation index and max contig prop (need N50, total length, maxlen)
	if [[ -n "$hmt_total_len" && -n "$hmt_n50" ]]; then
		hmt_frag_idx="$(awk -v TL="$hmt_total_len" -v N50="$hmt_n50" \
			'BEGIN{ if(TL>0){printf("%.4f", 1 - (N50/TL));} }')"
	fi
	if [[ -n "$hmt_total_len" && -n "$maxlen" ]]; then
		hmt_max_prop="$(awk -v TL="$hmt_total_len" -v ML="$maxlen" \
			'BEGIN{ if(TL>0){printf("%.4f", ML/TL);} }')"
	fi

	# CV of contig lengths from contig_table.tsv (Length column)
	if [[ -s "$ct" ]]; then
		# Column header assumed: Index ContigName ConservedGene GC Length Depth
		# We will read Length (5th) from data rows (NR>1).
		local awk_out
		awk_out="$(
			awk -F'\t' '
				NR==1{next}
				$1==""{next}
				{
					len[count]=$5+0; sum+=$5; sumsq+=$5*$5; count++
				}
				END{
					if(count>0){
						mean = sum/count
						var  = (sumsq - (sum*sum)/count)/count
						if(var<0) var=0
						sd   = sqrt(var)
						if(mean>0){
							cv = sd/mean
							printf "%.4f\n", cv
						}
					}
				}
			' "$ct"
		)"
		hmt_cv="$awk_out"
	fi
}

# Append OATK file paths as attributes under organelle "oatk"
# coverage label c10/c20/c30 etc.
_hifi_append_oatk_files() {
	local sp="$1" tier="$2" inum="$3" cov="$4" base_dir="$5"

	local label="c${cov}"
	local sdir="${sp}/${tier}/${inum}/summary-oatk-${cov}.txt"
	local assess_dir="${sp}/${tier}/${inum}/assess/oatk-${cov}/report"

	# --- File paths (mtDNA / ptDNA) ---
	local mt_gfa="${base_dir}/oatk-${cov}.mito.gfa"
	local mt_fa="${base_dir}/oatk-${cov}.mito.ctg.fasta"
	local mt_png="${base_dir}/oatk-${cov}.mito.png"
	local pt_gfa="${base_dir}/oatk-${cov}.pltd.gfa"
	local pt_fa="${base_dir}/oatk-${cov}.pltd.ctg.fasta"
	local pt_png="${base_dir}/oatk-${cov}.pltd.png"

	for k in mt_gfa mt_fa mt_png pt_gfa pt_fa pt_png; do
		local val="${!k:-}"
		local key="${label}_${k#*_}"
		if [[ -s "$val" ]]; then
			append_fact "$sp" "$tier" "$inum" "oatk" "attr" "$key" "$val"
		fi
	done

	# --- Summary metrics (memory/time/disk) ---
	local mem_kb="" mem_gb="" time_hms="" time_hours="" disk_gb=""
	if [[ -s "$sdir" ]]; then
		mem_kb=$(grep -E "Net increase" "$sdir" | awk '{print $3}' | head -n1)
		mem_gb=$(grep -E "Net increase" "$sdir" | sed -n 's/.*(\(.*GB\)).*/\1/p' | head -n1)
		time_hms=$(grep -E "Elapsed time" "$sdir" | awk '{print $3}' | head -n1)
		time_hours=$(grep -E "Elapsed time" "$sdir" | sed -n 's/.*(\(.*h\)).*/\1/p' | head -n1)
		disk_gb=$(grep -E "Disk used" "$sdir" | awk '{print $3}' | head -n1)
	fi
	append_fact "$sp" "$tier" "$inum" "oatk" "attr" "${label}_mem_kb" "$mem_kb"
	append_fact "$sp" "$tier" "$inum" "oatk" "attr" "${label}_mem_gb" "$mem_gb"
	append_fact "$sp" "$tier" "$inum" "oatk" "attr" "${label}_time_hms" "$time_hms"
	append_fact "$sp" "$tier" "$inum" "oatk" "attr" "${label}_time_hours" "$time_hours"
	append_fact "$sp" "$tier" "$inum" "oatk" "attr" "${label}_disk_gb" "$disk_gb"

	# --- HiMT metrics (4C) ---
	local geneset_completeness_prop="" num_contigs="" total_length="" N50=""
	local fragmentation_index="" max_contig_prop="" contig_length_cv=""

	if [[ -s "${assess_dir}/prefix.gene_integrity_summary.tsv" ]]; then
		geneset_completeness_prop=$(awk -F'\t' '$1=="geneset_completeness_prop"{print $2}' \
			"${assess_dir}/prefix.gene_integrity_summary.tsv" | head -n1)
	fi
	if [[ -s "${assess_dir}/prefix.mito_basic_info.tsv" ]]; then
		num_contigs=$(awk -F'\t' '$1=="Total contig number"{print $2}' \
			"${assess_dir}/prefix.mito_basic_info.tsv" | head -n1)
		total_length=$(awk -F'\t' '$1=="Total length"{print $2}' \
			"${assess_dir}/prefix.mito_basic_info.tsv" | head -n1)
		N50=$(awk -F'\t' '$1=="N50"{print $2}' \
			"${assess_dir}/prefix.mito_basic_info.tsv" | head -n1)
	fi
	if [[ -s "${assess_dir}/prefix.4c_metrics.tsv" ]]; then
		fragmentation_index=$(awk -F'\t' '$2=="fragmentation_index"{print $3}' \
			"${assess_dir}/prefix.4c_metrics.tsv" | head -n1)
		max_contig_prop=$(awk -F'\t' '$2=="max_contig_prop"{print $3}' \
			"${assess_dir}/prefix.4c_metrics.tsv" | head -n1)
		contig_length_cv=$(awk -F'\t' '$2=="contig_length_cv"{print $3}' \
			"${assess_dir}/prefix.4c_metrics.tsv" | head -n1)
	fi

	append_fact "$sp" "$tier" "$inum" "oatk" "attr" "${label}_geneset_completeness_prop" "$geneset_completeness_prop"
	append_fact "$sp" "$tier" "$inum" "oatk" "attr" "${label}_num_contigs" "$num_contigs"
	append_fact "$sp" "$tier" "$inum" "oatk" "attr" "${label}_total_length" "$total_length"
	append_fact "$sp" "$tier" "$inum" "oatk" "attr" "${label}_N50" "$N50"
	append_fact "$sp" "$tier" "$inum" "oatk" "attr" "${label}_fragmentation_index" "$fragmentation_index"
	append_fact "$sp" "$tier" "$inum" "oatk" "attr" "${label}_max_contig_prop" "$max_contig_prop"
	append_fact "$sp" "$tier" "$inum" "oatk" "attr" "${label}_contig_length_cv" "$contig_length_cv"
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

		# ===================== HiFi assembly entries ============================
		local root="${sp}/${TIER}/${INUM}"

		# ---- OATK variants: summary-oatk-*.txt -> single organelle "oatk" -----
		# if compgen -G "${root}/summary-oatk-*.txt" >/dev/null 2>&1; then
		if compgen -G "${root}/summary-oatk-*.txt" >/dev/null 2>&1; then
			for sumf in "${root}"/summary-oatk-30.txt; do
				[[ -e "$sumf" ]] || continue
				local cov tag_suffix
				cov="$(basename "$sumf" | sed -e 's/^summary-oatk-//' -e 's/\.txt$//')"
				tag_suffix="c${cov}"
				local oatk_dir="${root}/oatk"

				# local mtg="${oatk_dir}/oatk-${cov}.mito.gfa"
				# local mtf="${oatk_dir}/oatk-${cov}.mito.ctg.fasta"
				# local mtp="${oatk_dir}/oatk-${cov}.mito.png"
				# local ptg="${oatk_dir}/oatk-${cov}.pltd.gfa"
				# local ptf="${oatk_dir}/oatk-${cov}.pltd.ctg.fasta"
				# local ptp="${oatk_dir}/oatk-${cov}.pltd.png"

				# ----- OATK assemblies (c10,c20,c30 etc.) -----
				local oatk_dir="${sp}/${TIER}/${INUM}/oatk"
				if [[ -d "$oatk_dir" ]]; then
					# find all summary-oatk-*.txt and extract coverage numbers
					local cov
					for summary in "$sp/${TIER}/${INUM}"/summary-oatk-??.txt; do
						[[ -e "$summary" ]] || continue
						# summary-oatk-30.txt → 30
						cov="${summary##*-}" # 30.txt
						cov="${cov%%.txt}"   # 30
						_hifi_append_oatk_files "$sp" "$TIER" "$INUM" "$cov" "$oatk_dir"
						# you can also parse & append oatk resource metrics here if desired (mem/time/disk)
					done
				fi

			done
		else
			# No OATK summaries; still ensure there is an oatk block with placeholders
			append_fact "$sp" "$TIER" "$INUM" "oatk" "attr" "present" ""
		fi

		# ---- TIPPo assembly (single) ----
		local tip_dir="${root}/tippo"
		if [[ -d "$tip_dir" ]]; then
			# mt/pt GFA, FASTA, PNG (best-effort heuristics)
			local tip_mt_g tip_mt_f tip_mt_p tip_pt_g tip_pt_f tip_pt_p
			tip_mt_g="$(first_match "${tip_dir}"/*mitochondrial*flye/assembly_graph.gfa)"
			tip_mt_p="$(first_match "${tip_dir}"/*mitochondrial*flye/assembly_graph.png)"
			tip_pt_g="$(first_match "${tip_dir}"/*chloroplast*flye/assembly_graph.gfa)"
			tip_pt_f="$(first_match "${tip_dir}"/*organelle.chloroplast.fasta)"
			tip_pt_p="$(first_match "${tip_dir}"/*chloroplast*.png)"

			append_fact "$sp" "$TIER" "$INUM" "tippo" "file" "mt_gfa" "$tip_mt_g"
			# append_fact "$sp" "$TIER" "$INUM" "tippo" "file" "mt_fasta" "$tip_mt_f"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "file" "mt_png" "$tip_mt_p"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "file" "pt_gfa" "$tip_pt_g"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "file" "pt_fasta" "$tip_pt_f"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "file" "pt_png" "$tip_pt_p"

			# TIPPo summary
			_extract_summary_metrics "${root}/summary-tippo.txt"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "mem_kb" "${sum_mem_kb}"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "mem_gb" "${sum_mem_gb}"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "time_hms" "${sum_time_hms}"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "time_hours" "${sum_time_hours}"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "disk_gb" "${sum_disk_gb}"

			# TIPPo HiMT metrics
			_extract_himt_metrics "$sp" "$TIER" "$INUM" "tippo"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "geneset_completeness_prop" "${hmt_geneset}"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "num_contigs" "${hmt_contigs}"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "total_length" "${hmt_total_len}"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "N50" "${hmt_n50}"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "fragmentation_index" "${hmt_frag_idx}"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "max_contig_prop" "${hmt_max_prop}"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "contig_length_cv" "${hmt_cv}"
		else
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "present" ""
		fi

		# ---- HiMT assembly (self) ----
		local himt_dir="${root}/himt"
		if [[ -d "$himt_dir" ]]; then
			local himt_pt_g himt_pt_p himt_mt_g himt_mt_p himt_pt_f
			himt_pt_g="${himt_dir}/himt_chloroplast.gfa"
			himt_pt_p="${himt_dir}/himt_chloroplast.png"
			himt_mt_g="${himt_dir}/himt_mitochondrial.gfa"
			himt_mt_p="${himt_dir}/himt_mitochondrial.png"
			# choose path1 as primary pt fasta if exists
			himt_pt_f="$(first_match "${himt_dir}/chloroplast_path1.fa")"
			append_fact "$sp" "$TIER" "$INUM" "himt" "file" "mt_gfa" "$([[ -s "$himt_mt_g" ]] && echo "$himt_mt_g" || echo "")"
			append_fact "$sp" "$TIER" "$INUM" "himt" "file" "mt_fasta" ""
			append_fact "$sp" "$TIER" "$INUM" "himt" "file" "mt_png" "$([[ -s "$himt_mt_p" ]] && echo "$himt_mt_p" || echo "")"
			append_fact "$sp" "$TIER" "$INUM" "himt" "file" "pt_gfa" "$([[ -s "$himt_pt_g" ]] && echo "$himt_pt_g" || echo "")"
			append_fact "$sp" "$TIER" "$INUM" "himt" "file" "pt_fasta" "$himt_pt_f"
			append_fact "$sp" "$TIER" "$INUM" "himt" "file" "pt_png" "$([[ -s "$himt_pt_p" ]] && echo "$himt_pt_p" || echo "")"

			_extract_summary_metrics "${root}/summary-himt.txt"
			append_fact "$sp" "$TIER" "$INUM" "himt" "attr" "mem_kb" "${sum_mem_kb}"
			append_fact "$sp" "$TIER" "$INUM" "himt" "attr" "mem_gb" "${sum_mem_gb}"
			append_fact "$sp" "$TIER" "$INUM" "himt" "attr" "time_hms" "${sum_time_hms}"
			append_fact "$sp" "$TIER" "$INUM" "himt" "attr" "time_hours" "${sum_time_hours}"
			append_fact "$sp" "$TIER" "$INUM" "himt" "attr" "disk_gb" "${sum_disk_gb}"

			_extract_himt_metrics "$sp" "$TIER" "$INUM" "himt"
			append_fact "$sp" "$TIER" "$INUM" "himt" "attr" "geneset_completeness_prop" "${hmt_geneset}"
			append_fact "$sp" "$TIER" "$INUM" "himt" "attr" "num_contigs" "${hmt_contigs}"
			append_fact "$sp" "$TIER" "$INUM" "himt" "attr" "total_length" "${hmt_total_len}"
			append_fact "$sp" "$TIER" "$INUM" "himt" "attr" "N50" "${hmt_n50}"
			append_fact "$sp" "$TIER" "$INUM" "himt" "attr" "fragmentation_index" "${hmt_frag_idx}"
			append_fact "$sp" "$TIER" "$INUM" "himt" "attr" "max_contig_prop" "${hmt_max_prop}"
			append_fact "$sp" "$TIER" "$INUM" "himt" "attr" "contig_length_cv" "${hmt_cv}"
		else
			append_fact "$sp" "$TIER" "$INUM" "himt" "attr" "present" ""
		fi

		# ---- PMAT2 assembly ----
		local pmat_dir="${root}/pmat2/gfa_result"
		if [[ -d "$pmat_dir" ]]; then
			local pm_mt_g pm_mt_p pm_pt_g pm_pt_p pm_mt_f pm_pt_f
			pm_mt_g="${pmat_dir}/PMAT_mt_main.gfa"
			pm_mt_p="${pmat_dir}/PMAT_mt_main.png"
			pm_pt_g="${pmat_dir}/PMAT_pt_main.gfa"
			pm_pt_p="${pmat_dir}/PMAT_pt_main.png"
			pm_mt_f="${pmat_dir}/PMAT_mt.fa"
			pm_pt_f="${pmat_dir}/PMAT_pt.fa"

			append_fact "$sp" "$TIER" "$INUM" "pmat2" "file" "mt_gfa" "$([[ -s "$pm_mt_g" ]] && echo "$pm_mt_g" || echo "")"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "file" "mt_fasta" "$([[ -s "$pm_mt_f" ]] && echo "$pm_mt_f" || echo "")"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "file" "mt_png" "$([[ -s "$pm_mt_p" ]] && echo "$pm_mt_p" || echo "")"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "file" "pt_gfa" "$([[ -s "$pm_pt_g" ]] && echo "$pm_pt_g" || echo "")"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "file" "pt_fasta" "$([[ -s "$pm_pt_f" ]] && echo "$pm_pt_f" || echo "")"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "file" "pt_png" "$([[ -s "$pm_pt_p" ]] && echo "$pm_pt_p" || echo "")"

			_extract_summary_metrics "${root}/summary-pmat2.txt"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "attr" "mem_kb" "${sum_mem_kb}"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "attr" "mem_gb" "${sum_mem_gb}"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "attr" "time_hms" "${sum_time_hms}"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "attr" "time_hours" "${sum_time_hours}"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "attr" "disk_gb" "${sum_disk_gb}"

			_extract_himt_metrics "$sp" "$TIER" "$INUM" "pmat2"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "attr" "geneset_completeness_prop" "${hmt_geneset}"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "attr" "num_contigs" "${hmt_contigs}"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "attr" "total_length" "${hmt_total_len}"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "attr" "N50" "${hmt_n50}"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "attr" "fragmentation_index" "${hmt_frag_idx}"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "attr" "max_contig_prop" "${hmt_max_prop}"
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "attr" "contig_length_cv" "${hmt_cv}"
		else
			append_fact "$sp" "$TIER" "$INUM" "pmat2" "attr" "present" ""
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
	grep -E $'\t(data|short1data|short2data|oatk|tippo|himt|pmat2)\t' "$FACTS" || true
}

# assemble
python3 "${_POLAPLIB_DIR}/scripts/hifi_manifest_assemble.py" \
	--facts "$FACTS" \
	--set "$SET" --tier "$TIER" --inum "$INUM" \
	--out "$OUT" $([[ $PRETTY -eq 1 ]] && echo --pretty) \
	--species-codes "${SPECIES_CODES_TXT}"

cp "$FACTS" 1.tsv

[[ "$QUIET" -eq 1 ]] || echo "[INFO] Wrote: $OUT"
