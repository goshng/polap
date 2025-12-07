#!/usr/bin/env bash
#
# 2025-12-07
# polaplib/polap-bash-oatk-make-manifest.sh
#
# Version : v0.6.0  (2025-12-07)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Purpose:
#   Harvest per-species facts and assemble a structured JSON manifest via
#   scripts/oatk_manifest_assemble.py.
#
# This version (v0.6.0):
#   • Keeps legacy HiFi/polap-readassemble + OATK summary support.
#   • Adds harvesting of the new oatk–tiara–HiMT workflows laid out as:
#       <Species>/<tier>/<inum>/oatk-tiara-himt-<pattern>/{himt,tiara,oatk}/...
#     with corresponding logs at:
#       <Species>/<tier>/<inum>/stdout-oatk-tiara-himt-<pattern>-<step>.txt
#       <Species>/<tier>/<inum>/timing-oatk-tiara-himt-<pattern>-<step>.txt
#   • Parses `command time -v` timing logs for:
#       - Maximum resident set size (kbytes)
#       - Elapsed wall clock time (h:mm:ss or m:ss)
#       - User/System CPU times
#   • Emits per-pattern blocks under organelle names:
#       - organelle = "oatk-tiara-himt-<pattern>"
#     with keys like:
#       - himt_extract_fa, himt_reads_fq, himt_filter_mem_kb, ...
#       - tiara_out_tsv, tiara_mem_kb, ...
#       - c02_mito_gfa, c02_plastid_ctg_fa, c02_mem_kb, c02_time_hms, ...
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
# This code list should be from the manifest.
SPECIES_CODES_TXT="polap-species-codes-oatk.txt"
OUT="md/manifest.json"
INCLUDE_PTPT=0
PRETTY=0
DEBUG=0
QUIET="${QUIET:-1}"

usage() {
	cat <<EOF
Usage:
  $(basename "$0") --set some|auto|file:/path|Species_Name
                   --tier v2 --inum 0
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

first_match() {
	# Best-effort "find first file that matches"
	compgen -G "$1" >/dev/null 2>&1 && ls $1 2>/dev/null | head -n1 || true
}

# species list
get_species_list() {
	case "$SET" in
	auto)
		# auto-detect by presence of oatk-tiara-himt-* or legacy oatk directories
		{
			find . -mindepth 4 -maxdepth 4 -type d \
				-path "./*/${TIER}/${INUM}/oatk-tiara-himt-*"
			find . -mindepth 4 -maxdepth 4 -type d \
				-path "./*/${TIER}/${INUM}/oatk"
			find . -mindepth 4 -maxdepth 4 -type d \
				-path "./*/${TIER}/${INUM}/polap-readassemble"
		} 2>/dev/null |
			awk -F/ 'NF>=4{print $2}' |
			LC_ALL=C sort -u
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
	*)
		echo "$SET"
		;;
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

# Parse summary-*.txt (legacy) : set globals sum_mem_kb, sum_mem_gb, sum_time_hms, sum_time_hours, sum_disk_gb
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
	line="$(grep -m1 -E '^[[:space:]]*Disk[[:space:]]+used' "$f" || true)"
	sum_disk_gb=""
	if [[ -n "$line" ]]; then
		sum_disk_gb="$(
			printf '%s\n' "$line" |
				sed -E 's/.*Disk[[:space:]]+used:[^0-9]*([0-9]+(\.[0-9]+)?)[[:space:]]*GB.*/\1/'
		)"
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

# Parse "command time -v" stderr logs (timing-*.txt)
# Sets globals:
#   tv_mem_kb, tv_mem_gb, tv_time_hms, tv_time_hours, tv_user_s, tv_sys_s
_extract_timev_metrics() {
	local f="$1"
	tv_mem_kb=""
	tv_mem_gb=""
	tv_time_hms=""
	tv_time_hours=""
	tv_user_s=""
	tv_sys_s=""

	[[ ! -s "$f" ]] && return 0

	local out
	out="$(
		awk '
			/Maximum resident set size \(kbytes\)/{
				mem=$NF
			}
			/Elapsed \(wall clock\) time/{
				# Example:
				#   Elapsed (wall clock) time (h:mm:ss or m:ss): 0:31.54
				# Extract the token after the *last* colon ("0:31.54").
				if (match($0, /:[[:space:]]*([^[:space:]]+)$/, m)) {
					wall = m[1]
				}
			}
			/User time \(seconds\)/{
				user=$NF
			}
			/System time \(seconds\)/{
				sys=$NF
			}
			END{
				if (mem  == "") mem  = ""
				if (wall == "") wall = ""
				if (user == "") user = ""
				if (sys  == "") sys  = ""
				printf "%s\t%s\t%s\t%s\n", mem, wall, user, sys
			}
		' "$f"
	)"

	tv_mem_kb="$(printf '%s' "$out" | cut -f1)"
	tv_time_hms="$(printf '%s' "$out" | cut -f2)"
	tv_user_s="$(printf '%s' "$out" | cut -f3)"
	tv_sys_s="$(printf '%s' "$out" | cut -f4)"

	# Convert wall-clock to hours (handles h:mm:ss, m:ss, s, and fractional seconds)
	if [[ -n "$tv_time_hms" ]]; then
		tv_time_hours="$(
			awk -v t="$tv_time_hms" 'BEGIN{
				# Accept formats like:
				#   h:mm:ss      (e.g. 1:02:03)
				#   m:ss         (e.g. 0:31.54, 12:03)
				#   s or s.ss    (e.g. 31.54)
				n = split(t, a, ":");
				h = 0; m = 0; s = 0;
				if (n == 3) {
					h = a[1]; m = a[2]; s = a[3];
				} else if (n == 2) {
					m = a[1]; s = a[2];
				} else if (n == 1) {
					s = a[1];
				}
				sec = h*3600 + m*60 + s;
				if (sec > 0) {
					printf "%.3f", sec/3600.0;
				}
			}'
		)"
	else
		tv_time_hours=""
	fi

	# Convert kB to GB
	if [[ -n "$tv_mem_kb" ]]; then
		tv_mem_gb="$(
			awk -v kb="$tv_mem_kb" 'BEGIN{
				if (kb>0) { printf "%.3f", kb/1024.0/1024.0 } else { printf "" }
			}'
		)"
	else
		tv_mem_gb=""
	fi
}

# Append OATK file paths as attributes under organelle "oatk" (legacy summary-oatk-*.txt)
_oatk_append_oatk_files() {
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
# Harvest new oatk-tiara-himt-* layouts for a species
#   <sp>/<tier>/<inum>/oatk-tiara-himt-<pattern>/{himt,tiara,oatk}
# ------------------------------------------------------------------------------
harvest_oatk_tiara_himt() {
	local sp="$1" tier="$2" inum="$3"
	local root="${sp}/${tier}/${inum}"

	shopt -s nullglob
	local od
	for od in "${root}"/oatk-tiara-himt-*; do
		[[ -d "$od" ]] || continue
		local b pattern org
		b="$(basename "$od")"           # e.g. oatk-tiara-himt-x-h
		pattern="${b#oatk-tiara-himt-}" # e.g. x-h
		org="oatk-tiara-himt-${pattern}"

		# Make sure the block exists
		append_fact "$sp" "$tier" "$inum" "$org" "attr" "present" "1"

		# ---------------------- HiMT (filter) ----------------------
		local himt_dir="${od}/himt"
		if [[ -d "$himt_dir" ]]; then
			local h_extract="${himt_dir}/extract.fa"
			local h_fq="${himt_dir}/himt.fq.gz"
			local h_proc="${himt_dir}/process.fa"
			local h_prop="${himt_dir}/proportion.fa"

			[[ -s "$h_extract" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "himt_extract_fa" "$h_extract"
			[[ -s "$h_fq" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "himt_reads_fq" "$h_fq"
			[[ -s "$h_proc" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "himt_process_fa" "$h_proc"
			[[ -s "$h_prop" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "himt_proportion_fa" "$h_prop"

			local hf_stdout="${root}/stdout-oatk-tiara-himt-${pattern}-himt-filter.txt"
			local hf_timing="${root}/timing-oatk-tiara-himt-${pattern}-himt-filter.txt"
			[[ -s "$hf_stdout" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "himt_filter_stdout" "$hf_stdout"
			[[ -s "$hf_timing" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "himt_filter_timing" "$hf_timing"

			if [[ -s "$hf_timing" ]]; then
				_extract_timev_metrics "$hf_timing"
				append_fact "$sp" "$tier" "$inum" "$org" "attr" "himt_filter_mem_kb" "$tv_mem_kb"
				append_fact "$sp" "$tier" "$inum" "$org" "attr" "himt_filter_mem_gb" "$tv_mem_gb"
				append_fact "$sp" "$tier" "$inum" "$org" "attr" "himt_filter_time_hms" "$tv_time_hms"
				append_fact "$sp" "$tier" "$inum" "$org" "attr" "himt_filter_time_hours" "$tv_time_hours"
				append_fact "$sp" "$tier" "$inum" "$org" "attr" "himt_filter_user_s" "$tv_user_s"
				append_fact "$sp" "$tier" "$inum" "$org" "attr" "himt_filter_sys_s" "$tv_sys_s"
			fi
		fi

		# ---------------------- Tiara ------------------------------
		local tiara_dir="${od}/tiara"
		if [[ -d "$tiara_dir" ]]; then
			local t_out="${tiara_dir}/tiara.out.tsv"
			local t_mito_in="${tiara_dir}/mitochondrion_tiara.input.fa"
			local t_pltd_in="${tiara_dir}/plastid_tiara.input.fa"
			local t_mito_fq="${tiara_dir}/tiara.mito.fq.gz"
			local t_pltd_fq="${tiara_dir}/tiara.plastid.fq.gz"

			[[ -s "$t_out" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "tiara_out_tsv" "$t_out"
			[[ -s "$t_mito_in" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "tiara_mito_input_fa" "$t_mito_in"
			[[ -s "$t_pltd_in" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "tiara_plastid_input_fa" "$t_pltd_in"
			[[ -s "$t_mito_fq" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "tiara_mito_fq" "$t_mito_fq"
			[[ -s "$t_pltd_fq" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "tiara_plastid_fq" "$t_pltd_fq"

			local t_stdout="${root}/stdout-oatk-tiara-himt-${pattern}-tiara.txt"
			local t_timing="${root}/timing-oatk-tiara-himt-${pattern}-tiara.txt"
			[[ -s "$t_stdout" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "tiara_stdout" "$t_stdout"
			[[ -s "$t_timing" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "tiara_timing" "$t_timing"

			if [[ -s "$t_timing" ]]; then
				_extract_timev_metrics "$t_timing"
				append_fact "$sp" "$tier" "$inum" "$org" "attr" "tiara_mem_kb" "$tv_mem_kb"
				append_fact "$sp" "$tier" "$inum" "$org" "attr" "tiara_mem_gb" "$tv_mem_gb"
				append_fact "$sp" "$tier" "$inum" "$org" "attr" "tiara_time_hms" "$tv_time_hms"
				append_fact "$sp" "$tier" "$inum" "$org" "attr" "tiara_time_hours" "$tv_time_hours"
				append_fact "$sp" "$tier" "$inum" "$org" "attr" "tiara_user_s" "$tv_user_s"
				append_fact "$sp" "$tier" "$inum" "$org" "attr" "tiara_sys_s" "$tv_sys_s"
			fi
		fi

		# ---------------------- OATK (per coverage) ----------------
		local oatk_dir="${od}/oatk"
		if [[ -d "$oatk_dir" ]]; then
			local mt_gfa
			for mt_gfa in "${oatk_dir}"/oatk-"${pattern}"-*.mito.gfa; do
				[[ -e "$mt_gfa" ]] || continue
				local base cov label
				base="$(basename "$mt_gfa")"         # oatk-x-h-02.mito.gfa
				local tmp="${base#oatk-${pattern}-}" # 02.mito.gfa
				cov="${tmp%%.mito.gfa}"              # 02
				label="c${cov}"

				local mt_png="${oatk_dir}/oatk-${pattern}-${cov}.mito.png"
				local mt_ctg="${oatk_dir}/oatk-${pattern}-${cov}.mito.ctg.fasta"
				local pt_gfa="${oatk_dir}/oatk-${pattern}-${cov}.pltd.gfa"
				local pt_png="${oatk_dir}/oatk-${pattern}-${cov}.pltd.png"
				local pt_ctg="${oatk_dir}/oatk-${pattern}-${cov}.pltd.ctg.fasta"

				[[ -s "$mt_gfa" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "${label}_mito_gfa" "$mt_gfa"
				[[ -s "$mt_png" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "${label}_mito_png" "$mt_png"
				[[ -s "$mt_ctg" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "${label}_mito_ctg_fa" "$mt_ctg"
				[[ -s "$pt_gfa" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "${label}_plastid_gfa" "$pt_gfa"
				[[ -s "$pt_png" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "${label}_plastid_png" "$pt_png"
				[[ -s "$pt_ctg" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "${label}_plastid_ctg_fa" "$pt_ctg"

				local o_stdout="${root}/stdout-oatk-tiara-himt-${pattern}-c${cov}.txt"
				local o_timing="${root}/timing-oatk-tiara-himt-${pattern}-c${cov}.txt"
				[[ -s "$o_stdout" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "${label}_stdout" "$o_stdout"
				[[ -s "$o_timing" ]] && append_fact "$sp" "$tier" "$inum" "$org" "file" "${label}_timing" "$o_timing"

				if [[ -s "$o_timing" ]]; then
					_extract_timev_metrics "$o_timing"
					append_fact "$sp" "$tier" "$inum" "$org" "attr" "${label}_mem_kb" "$tv_mem_kb"
					append_fact "$sp" "$tier" "$inum" "$org" "attr" "${label}_mem_gb" "$tv_mem_gb"
					append_fact "$sp" "$tier" "$inum" "$org" "attr" "${label}_time_hms" "$tv_time_hms"
					append_fact "$sp" "$tier" "$inum" "$org" "attr" "${label}_time_hours" "$tv_time_hours"
					append_fact "$sp" "$tier" "$inum" "$org" "attr" "${label}_user_s" "$tv_user_s"
					append_fact "$sp" "$tier" "$inum" "$org" "attr" "${label}_sys_s" "$tv_sys_s"
				fi
			done
		fi
	done
	shopt -u nullglob || true
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
			local fq
			fq="$(first_match "${sdir}/l.fq.gz")"
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
				local sfq
				sfq="$(first_match "${sdir}/s_${mate}.fq.gz")"
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

		# ------------------------------ PT files (legacy polap-readassemble) ----
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

		# ------------------------------ MT files (legacy polap-readassemble) ----
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

		extract_pattern_cov() {
			local f="$1"
			local base mid pattern cov
			base="${f##*/}"
			mid="${base#timing-oatk-tiara-himt-}"
			mid="${mid%.txt}"

			pattern="${mid%-c*}" # e.g. h-t
			cov="${mid##*-c}"    # e.g. 02

			printf "%s-%s\n" "$pattern" "$cov"
		}

		extract_pattern_only() {
			local f="$1"
			local base mid pattern
			base="${f##*/}"
			mid="${base#timing-oatk-tiara-himt-}"
			mid="${mid%.txt}"
			pattern="${mid%-c*}"
			printf "%s\n" "$pattern"
		}

		# ===================== legacy oatk assembly entries =====================
		local root="${sp}/${TIER}/${INUM}"

		if compgen -G "${root}/timing-oatk-tiara-himt-*.txt" >/dev/null 2>&1; then
			local summary
			for summary in "$root"/timing-oatk-tiara-himt-*.txt; do
				[[ -e "$summary" ]] || continue
				local cov=$(extract_pattern_cov "$summary")
				local pattern=$(extract_pattern_only "$summary")
				echo "summary: $summary" >&2
				echo "cov: $cov" >&2
				echo "pattern: $pattern" >&2
				local oatk_dir="${root}/oatk-tiara-himt-${pattern}/oatk"
				if [[ -d "$oatk_dir" ]]; then
					_oatk_append_oatk_files "$sp" "$TIER" "$INUM" "$cov" "$oatk_dir"
				fi
			done
		else
			# No legacy OATK summaries
			append_fact "$sp" "$TIER" "$INUM" "oatk" "attr" "present" ""
		fi

		# ---- TIPPo assembly (legacy) ----
		local tip_dir="${root}/tippo"
		if [[ -d "$tip_dir" ]]; then
			local tip_mt_g tip_mt_p tip_pt_g tip_pt_f tip_pt_p
			tip_mt_g="$(first_match "${tip_dir}"/*mitochondrial*flye/assembly_graph.gfa)"
			tip_mt_p="$(first_match "${tip_dir}"/*mitochondrial*flye/assembly_graph.png)"
			tip_pt_g="$(first_match "${tip_dir}"/*chloroplast*flye/assembly_graph.gfa)"
			tip_pt_f="$(first_match "${tip_dir}"/*organelle.chloroplast.fasta)"
			tip_pt_p="$(first_match "${tip_dir}"/*chloroplast*.png)"

			append_fact "$sp" "$TIER" "$INUM" "tippo" "file" "mt_gfa" "$tip_mt_g"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "file" "mt_png" "$tip_mt_p"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "file" "pt_gfa" "$tip_pt_g"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "file" "pt_fasta" "$tip_pt_f"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "file" "pt_png" "$tip_pt_p"

			_extract_summary_metrics "${root}/summary-tippo.txt"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "mem_kb" "${sum_mem_kb}"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "mem_gb" "${sum_mem_gb}"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "time_hms" "${sum_time_hms}"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "time_hours" "${sum_time_hours}"
			append_fact "$sp" "$TIER" "$INUM" "tippo" "attr" "disk_gb" "${sum_disk_gb}"

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

		# ---- HiMT assembly (legacy) ----
		local himt_dir="${root}/himt"
		if [[ -d "$himt_dir" ]]; then
			local himt_pt_g himt_pt_p himt_mt_g himt_mt_p himt_pt_f
			himt_pt_g="${himt_dir}/himt_chloroplast.gfa"
			himt_pt_p="${himt_dir}/himt_chloroplast.png"
			himt_mt_g="${himt_dir}/himt_mitochondrial.gfa"
			himt_mt_p="${himt_dir}/himt_mitochondrial.png"
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

		# ---- PMAT2 assembly (legacy) ----
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

		# ===================== NEW oatk-tiara-himt patterns =====================
		harvest_oatk_tiara_himt "$sp" "$TIER" "$INUM"

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
	grep -E $'\t(data|short1data|short2data|oatk|tippo|himt|pmat2|oatk-tiara-himt-)\t' "$FACTS" || true
}

# assemble
python3 "${_POLAPLIB_DIR}/scripts/oatk_manifest_assemble.py" \
	--facts "$FACTS" \
	--set "$SET" --tier "$TIER" --inum "$INUM" \
	--out "$OUT" $([[ $PRETTY -eq 1 ]] && echo --pretty) \
	--species-codes "${SPECIES_CODES_TXT}"

cp "$FACTS" 1.tsv

[[ "$QUIET" -eq 1 ]] || echo "[INFO] Wrote: $OUT"
