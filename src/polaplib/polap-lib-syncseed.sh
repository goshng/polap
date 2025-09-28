################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

################################################################################
# Convert numbers between different units.
################################################################################

################################################################################
# Ensure that the current script is sourced only once
source "${_POLAPLIB_DIR}/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u
if [[ -n "${!_POLAP_INCLUDE_}" ]]; then
	set -u
	return 0
fi
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
	echo "[ERROR] This script must be sourced, not executed: use 'source $BASH_SOURCE'" >&2
	return 1 2>/dev/null || exit 1
fi
: "${_POLAP_DEBUG:=0}"
: "${_POLAP_RELEASE:=0}"

get_top_syncasm_opts() {
	local root="${1:-run_auto}"
	local verbose="${2:-0}" # 0: quiet (default), nonzero: print top-10 table

	# Ensure global associative array
	declare -gA SYNCASM_TOP_OPTS
	SYNCASM_TOP_OPTS=()

	# Pattern: exp_aN_kN_sN_cN_wxN_uN (numbers of any length)
	local pat='^exp_a[0-9]+_k[0-9]+_s[0-9]+_c[0-9]+_wx[0-9]+_u[0-9]+$'

	# Build table: name utgs total_bp
	local table
	table=$(
		shopt -s nullglob
		for d in "$root"/*/; do
			local nm="$(basename "$d")"
			[[ "$nm" =~ $pat ]] || continue

			local gfa="${d%/}/syncasm.asm.utg.final.gfa"
			[[ -s "$gfa" ]] || continue

			local n L
			n=$(awk '$1=="S"{c++} END{print c+0}' "$gfa")
			L=$(awk '
        $1=="S"{
          if($3!="*") len+=length($3)
          else for(i=4;i<=NF;i++) if($i~/^LN:i:/){split($i,a,":");len+=a[3]}
        }
        END{print len+0}
      ' "$gfa")

			printf "%s\t%s\t%s\n" "$nm" "$n" "$L"
		done | sort -k3,3nr -k2,2n
	) || return 1

	[[ -n "$table" ]] || {
		echo "# get_top_syncasm_opts: no matching runs under '$root'"
		return 2
	}

	# Best row
	local top name
	top=$(head -n1 <<<"$table")
	name=$(awk '{print $1}' <<<"$top")

	# Optional top-10 table
	if ((verbose)); then
		echo "# Top 10 runs (sorted by total_bp desc, utgs asc):"
		printf "%-3s  %-40s  %-8s  %-12s\n" "No" "name" "utgs" "total_bp"
		local i=0
		while IFS=$'\t' read -r nm utgs bp; do
			((i++))
			printf "%-3d  %-40s  %-8s  %-12s\n" "$i" "$nm" "$utgs" "$bp"
			((i >= 10)) && break
		done <<<"$table"
		echo
	fi

	# Parse options into SYNCASM_TOP_OPTS
	if [[ "$name" =~ _k([0-9]+) ]]; then SYNCASM_TOP_OPTS[k]="${BASH_REMATCH[1]}"; fi
	if [[ "$name" =~ _s([0-9]+) ]]; then SYNCASM_TOP_OPTS[s]="${BASH_REMATCH[1]}"; fi
	if [[ "$name" =~ _c([0-9]+) ]]; then SYNCASM_TOP_OPTS[c]="${BASH_REMATCH[1]}"; fi
	if [[ "$name" =~ _u([0-9]+) ]]; then SYNCASM_TOP_OPTS[u]="${BASH_REMATCH[1]}"; fi
	if [[ "$name" =~ _a([0-9]+) ]]; then
		SYNCASM_TOP_OPTS[a]="0.${BASH_REMATCH[1]}"
	fi
	if [[ "$name" =~ _wx([0-9]+) ]]; then
		SYNCASM_TOP_OPTS[wx]="0.${BASH_REMATCH[1]}"
	fi

	# Report chosen best
	echo "# Top run: $name"
	for key in "${!SYNCASM_TOP_OPTS[@]}"; do
		printf "  %s = %s\n" "$key" "${SYNCASM_TOP_OPTS[$key]}"
	done
}

recreate_syncasm_name() {
	local a="${SYNCASM_TOP_OPTS[a]#0.}"
	a=$((10#$a))
	a=$(printf "%03d" "$a")

	local wx="${SYNCASM_TOP_OPTS[wx]#0.}"
	wx=$((10#$wx))
	wx=$(printf "%03d" "$wx")

	local k="${SYNCASM_TOP_OPTS[k]}"
	local s="${SYNCASM_TOP_OPTS[s]}"
	local c="${SYNCASM_TOP_OPTS[c]}"
	local u="${SYNCASM_TOP_OPTS[u]}"

	echo "exp_a${a}_k${k}_s${s}_c${c}_wx${wx}_u${u}"
}

# Keep ONT reads that do NOT map to the reference (based on PAF)
# Usage:
#   filter_ont_unmapped_by_ref ref.fa reads.fastq[.gz] outprefix \
#       [min_identity] [min_qcov] [min_tcov] [min_mapq] [threads]
#
# Defaults:
#   min_identity = 0.85        # col10/col11  (matches / aln_len)
#   min_qcov     = 0.50        # (qend - qstart) / qlen
#   min_tcov     = 0.00        # (tend - tstart) / tlen
#   min_mapq     = 0           # PAF col12
#   threads      = 8
#
# Outputs:
#   outprefix.paf              : PAF alignments
#   outprefix.mapped.ids       : read IDs that pass the filters
#   outprefix.unmapped.ids     : read IDs NOT passing the filters (complement)
#   outprefix.unmapped.fastq.gz: unmapped reads
filter_ont_unmapped_by_ref() {
	local ref="$1"
	local reads="$2"
	local out="$3"
	# ONT: 0.65
	local min_id="${4:-0.65}"
	local min_qcov="${5:-0.50}"
	local min_tcov="${6:-0.00}"
	local min_mapq="${7:-0}"
	local threads="${8:-8}"

	if [[ -z "$ref" || -z "$reads" || -z "$out" ]]; then
		echo "Usage: filter_ont_unmapped_by_ref ref.fa reads.fq[.gz] outprefix [id] [qcov] [tcov] [mapq] [threads]" >&2
		return 2
	fi

	# 1) Map reads → PAF (no secondary alignments)

	minimap2 -t "$threads" -x "${_arg_minimap2_data_type}" --secondary=no "$ref" "$reads" >"${out}.paf"

	# 2) Collect mapped IDs that PASS thresholds
	#    PAF columns: 1=qname, 2=qlen, 3=qstart, 4=qend, 6=tname, 7=tlen, 8=tstart, 9=tend, 10=nmatch, 11=aln_len, 12=mapq
	awk -v ID="$min_id" -v QC="$min_qcov" -v TC="$min_tcov" -v MQ="$min_mapq" '
    $10>0 && $11>0 {
      id   = $10/$11
      qcov = ($4-$3)/$2
      tcov = ($9-$8)/$7
      mapq = $12+0
      if (id>=ID && qcov>=QC && tcov>=TC && mapq>=MQ) print $1
    }
  ' "${out}.paf" | sort -u >"${out}.mapped.ids"

	# 3) Collect ALL read IDs from the FASTQ (gz or plain)
	if [[ "$reads" =~ \.gz$ ]]; then
		zcat -- "$reads" |
			awk 'NR%4==1{split($1,a," "); print substr(a[1],2)}' |
			sort -u >"${out}.all.ids"
	else
		awk 'NR%4==1{split($1,a," "); print substr(a[1],2)}' "$reads" |
			sort -u >"${out}.all.ids"
	fi
	# 4) Unmapped = ALL minus MAPPED
	comm -23 "${out}.all.ids" "${out}.mapped.ids" >"${out}.unmapped.ids"

	# 5) Extract unmapped reads
	seqtk subseq "$reads" "${out}.unmapped.ids" | gzip -c >"${out}.unmapped.fastq.gz"
	# seqtk subseq "$reads" "${out}.unmapped.ids" >"${out}.unmapped.fastq"

	# 6) Tiny report
	echo "# mapped ids   : $(wc -l <"${out}.mapped.ids")"
	echo "# all ids      : $(wc -l <"${out}.all.ids")"
	echo "# unmapped ids : $(wc -l <"${out}.unmapped.ids")"
	echo "# output       : ${out}.unmapped.fastq.gz"
}

_polap_lib_syncseed-run() {
	local s="${1:-1}"
	_polap_log0 "Make it!"

	#
	# preprocessing
	#
	local READS="${_arg_long_reads}"
	local OUT="${_arg_outdir}"

	local syncasmdir="${_arg_outdir}/syncasm"
	# remove ptDNA reads
	if ((s == 0)); then
		_polap_lib_conda-ensure_conda_env polap || exit 1
		_polap_log0 "remove ptDNA reads: filter out ptDNA reads by mapping"
		rm -rf "${syncasmdir}/s0"
		mkdir -p "${syncasmdir}/s0"
		filter_ont_unmapped_by_ref "${_arg_unpolished_fasta}" \
			"${_arg_long_reads}" \
			"${syncasmdir}/s0/pt"
		conda deactivate
	fi

	if ((s == 1)); then
		_polap_lib_conda-ensure_conda_env polap-ont || exit 1
		bash "${_POLAPLIB_DIR}/polap-bash-oatk-ont-s1.sh" \
			--trim \
			--scrub \
			--lenfilt \
			--reads "$READS" \
			--out "$OUT" \
			--oatkdb "$HOME/OatkDB"
		conda deactivate
	fi

	# check -c from top to 1.
	# see if we find the longest mtDNA seed contigs
	# minimap2
	# seqtk subseq
	# flye
	# polap-oga
	if ((s == 2)); then
		_polap_log0 "select mtDNA seed contigs"

		_polap_lib_conda-ensure_conda_env polap-oatk || exit 1

		# ONT default tech; no HPC scaling; keep your k/s
		rm -rf "${syncasmdir}/s2"
		mkdir -p "${syncasmdir}/s2"
		python "${_POLAPLIB_DIR}/polap-py-syncasm-auto.py" \
			--out "${syncasmdir}/s2" \
			--reads "${syncasmdir}/s0/pt.unmapped.fastq.gz" \
			--tech ont \
			--hpc \
			--k 91 --s 21 \
			-t 16 >"${syncasmdir}/s2/summary.txt"

		get_top_syncasm_opts "${syncasmdir}/s2" 1

		# Access keys like an associative array
		_polap_log0 "Round 1"
		_polap_log0 "Best k = ${SYNCASM_TOP_OPTS[k]}"
		_polap_log0 "Best s = ${SYNCASM_TOP_OPTS[s]}"
		_polap_log0 "Best c = ${SYNCASM_TOP_OPTS[c]}"
		_polap_log0 "Best a = ${SYNCASM_TOP_OPTS[a]}"
		_polap_log0 "Best weak-cross = ${SYNCASM_TOP_OPTS[wx]}"
		_polap_log0 "Best unzip-round = ${SYNCASM_TOP_OPTS[u]}"

		local c

		mkdir -p "${syncasmdir}/s3"
		# for c in 100 50 30 20 10 5 4 3 2 1; do
		for c in 10 5 4 3 2 1; do
			SYNCASM_TOP_OPTS[c]="$c"
			outname=$(recreate_syncasm_name)
			_polap_log0 "starting $outname"

			mkdir -p "${syncasmdir}/s3/$outname"
			syncasm -c "$c" \
				-k "${SYNCASM_TOP_OPTS[k]}" \
				-s "${SYNCASM_TOP_OPTS[s]}" \
				-a "${SYNCASM_TOP_OPTS[a]}" \
				--weak-cross "${SYNCASM_TOP_OPTS[wx]}" \
				--unzip-round "${SYNCASM_TOP_OPTS[u]}" \
				-t "${_arg_threads}" \
				-o "${syncasmdir}/s3/$outname/syncasm.asm" \
				"${syncasmdir}/s2/reads.hpc.fastq" \
				>"${syncasmdir}/s3/$outname/job.out" \
				2>"${syncasmdir}/s3/$outname/job.err"
		done

		_polap_log0 "All syncasm jobs finished."
		get_top_syncasm_opts "${syncasmdir}/s3" 1

		# Access keys like an associative array
		_polap_log0 "Round 2"
		_polap_log0 "Best k = ${SYNCASM_TOP_OPTS[k]}"
		_polap_log0 "Best s = ${SYNCASM_TOP_OPTS[s]}"
		_polap_log0 "Best c = ${SYNCASM_TOP_OPTS[c]}"
		_polap_log0 "Best a = ${SYNCASM_TOP_OPTS[a]}"
		_polap_log0 "Best weak-cross = ${SYNCASM_TOP_OPTS[wx]}"
		_polap_log0 "Best unzip-round = ${SYNCASM_TOP_OPTS[u]}"

		#
		outname=$(recreate_syncasm_name)
		# "${syncasmdir}/s3/$outname/syncasm.asm.utg.final.gfa"

		conda deactivate
	fi

	if ((s == 3)); then
		_polap_lib_conda-ensure_conda_env polap || exit 1
		# annotate to select mt seed contigs

		conda deactivate
	fi

	if ((s == 4)); then
		_polap_lib_conda-ensure_conda_env polap-oatk || exit 1
		bash "${_POLAPLIB_DIR}/polap-bash-oatk-ont-s2.sh" \
			--reads-pre "$READS" \
			--out "$OUT" \
			--oatkdb "$HOME/OatkDB" \
			--k-list 121 \
			--smer 31 \
			--a 0.25 \
			--c 10 \
			--max-bubble 200000 \
			--max-tip 20000 \
			--weak-cross 0.30 \
			--unzip-round 2 \
			--threads 32
		conda deactivate
	fi

	return
	# bbmap
	# meryl
	# seqtk
	# seqkit

	#############################################
	# Stage-1: assemble (backbone only)
	#############################################
	log 1 "[stage1] assemble backbone (k=${K_LIST_A}; smer=${SMER}; -a=${AARC})"
	bash "$ONT_SCRIPT" \
		--reads "$CUR_ASM" \
		--out "$OUT/ont_stage1" \
		--oatkdb "$OATKDB" --clade "$CLADE" \
		--threads "$THREADS" \
		"${HPC_FLAG[@]}" \
		--k-list "$K_LIST_A" --smer "$SMER" --a "$AARC" --weak-cross "$WEAKX" --unzip-round "$UNZIP" \
		$([[ $NO_EC -eq 1 ]] && echo --no-read-ec) \
		-s qc,assemble,summary -v

	BACKBONE="$OUT/ont_stage1/k1/unitigs.fa"
	[[ -s "$BACKBONE" ]] || die "backbone unitigs not found: $BACKBONE"

	#############################################
	# Bait: backbone-based k-mer bait → subset
	#############################################
	mkdir -p bait
	SUBSET="bait/ont.mt.fq"

	if [[ "$BAIT_METHOD" == "bbduk" ]]; then
		need "$BBDUK_BIN"
		log 1 "[bait] bbduk: ref=$BACKBONE → $SUBSET"
		"$BBDUK_BIN" in="$CUR" outm="$SUBSET" outu=bait/nonmt.fq ref="$BACKBONE" k=31 hdist=1 threads="$THREADS"
	else
		need meryl
		need meryl-lookup
		log 1 "[bait] meryl k=$MERYL_K"
		meryl k="$MERYL_K" count "$BACKBONE" output bait/seeds.k"$MERYL_K".meryl
		meryl-lookup -existence -sequence "$CUR" bait/seeds.k"$MERYL_K".meryl |
			awk '/^>/{print substr($0,2)}' >bait/mt.ids
		seqkit grep -f bait/mt.ids "$CUR" -o "$SUBSET"
	fi

	log 1 "[bait] subset size: $( (wc -c <"$SUBSET") 2>/dev/null || echo 0) bytes"

	#############################################
	# Stage-2: lift OR polish → annotate → PF → summary
	#############################################
	if [[ $USE_LIFT -eq 1 ]]; then
		# Lifter (RLE) + 1× polish after lift
		log 1 "[stage2] LIFT (RLE) + 1× polish (fast), then annotate & PF"
		bash "$ONT_SCRIPT" \
			--reads "$CUR_ASM" \
			--reads-map "$SUBSET" \
			--map-chunks "$MAP_CHUNKS" --mm2-extra "$MM2_EXTRA" \
			--out "$OUT/ont_stage2" \
			--oatkdb "$OATKDB" --clade "$CLADE" \
			--threads "$THREADS" \
			--lift --polish-after-lift \
			--racon-rounds 1 \
			$([[ $NO_MEDAKA -eq 1 ]] && echo --no-medaka || echo --medaka-model "$MEDAKA_MODEL") \
			-s lift,polish,annotate,pathfinder,summary -v
	else
		# Classic polish on subset
		log 1 "[stage2] POLISH (racon/medaka) on subset, then annotate & PF"
		bash "$ONT_SCRIPT" \
			--reads "$CUR_ASM" \
			--reads-map "$SUBSET" \
			--map-chunks "$MAP_CHUNKS" --mm2-extra "$MM2_EXTRA" \
			--out "$OUT/ont_stage2" \
			--oatkdb "$OATKDB" --clade "$CLADE" \
			--threads "$THREADS" \
			--racon-rounds "$RACON_ROUNDS" \
			$([[ $NO_MEDAKA -eq 1 ]] && echo --no-medaka || echo --medaka-model "$MEDAKA_MODEL") \
			-s polish,annotate,pathfinder,summary -v
	fi

	#############################################
	# Side-report
	#############################################
	log 1 "[EMIT] sidekicks report"
	{
		echo -e "key\tvalue"
		echo -e "reads_in\t$READS"
		echo -e "out\t$OUT"
		echo -e "hpc\t$([[ $HPC_ENABLE -eq 1 ]] && echo ON || echo OFF)"
		echo -e "trim\t$([[ $DO_TRIM -eq 1 ]] && echo ON || echo OFF)"
		echo -e "scrub\t$([[ $DO_SCRUB -eq 1 ]] && echo ON || echo OFF)"
		echo -e "len_min\t$([[ $DO_LENFILT -eq 1 ]] && echo $LEN_MIN || echo .)"
		echo -e "keep_pct\t$([[ $DO_LENFILT -eq 1 && $(
			command -v filtlong >/dev/null
			echo $?
		) -eq 0 ]] && echo $KEEP_PCT || echo .)"
		echo -e "mt_bait_prefilter\t$([[ $DO_MT_BAIT -eq 1 ]] && echo ON || echo OFF)"
		echo -e "drop_plastid_prefilter\t$([[ $DO_DROP_PLASTID -eq 1 ]] && echo ON || echo OFF)"
		echo -e "duplex_only\t$([[ $DO_DUPLEX_ONLY -eq 1 ]] && echo $DUPLEX_TAG || echo OFF)"
		echo -e "rare_mask\t$([[ $DO_RARE_MASK -eq 1 ]] && echo K=$RARE_K,max=$RARE_MAX || echo OFF)"
		echo -e "corrector\t$([[ $DO_CORRECT -eq 1 ]] && echo CANU || ([[ $DO_CORRECT -eq 2 ]] && echo RATATOSK || echo OFF))"
		echo -e "bait_method\t$BAIT_METHOD"
		echo -e "map_chunks\t$MAP_CHUNKS"
		echo -e "mm2_extra\t${MM2_EXTRA:-.}"
		echo -e "ont_script\t$ONT_SCRIPT"
		echo -e "k_list_stage1\t$K_LIST_A"
		echo -e "racon_rounds\t$RACON_ROUNDS"
		echo -e "medaka_model\t$([[ $NO_MEDAKA -eq 1 ]] && echo OFF || echo $MEDAKA_MODEL)"
	} >"$OUT/sidekicks.report.tsv"

	log 1 "[done] Two-stage ONT run completed."
	log 1 "  Stage-1 backbone : $OUT/ont_stage1/k1/unitigs.fa"
	log 1 "  Stage-2 outputs  : $OUT/ont_stage2 (polished/lifted, annotated, pathfinder, summary)"

	return

	# seqtk hpc
	local CUR="${OUT}"
	local HPC_FILE="${HPC_OUT:-$OUT/reads.hpc.fa}"
	_polap_log1 "[hpc] seqtk hpc → $HPC_FILE"
	seqtk hpc "$CUR" >"$HPC_FILE"

	# syncasm
	oatk -k 1001 -s 31 -o $OUT/oatk.asm

	# minimap2
	minimap2

	# seqtk subseq
	# flye
	# polap-oga

	conda deactivate
}
