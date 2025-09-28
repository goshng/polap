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
# Functions for subcommand template ...
# Describe what they are and what they do.
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

function _run_polap_syncassemble {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	local polap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  polap ${polap_cmd} - filter out nuclear-derived PacBio HiFi reads

Synopsis:
  polap ${polap_cmd} [options]

Description:
  polap ${polap_cmd} uses Oatk's syncasm to filter out nuclear-derived reads.

Inputs:
  PacBio HiFi long-read

Outputs:
  Plastid genome assembly and DNA sequence:
    ${_arg_outdir}/pt-sync.0.gfa
    ${_arg_outdir}/pt-sync.0.fa

  Mitochondrial genome assembly:
    ${_arg_outdir}/mt-sync.0.gfa
  
Options:
  -l FASTQ: reads data file

Examples:
  Assemble plant organelle genomes from PacBio HiFi long-read data:
    polap ${polap_cmd} --pacbio-hifi -l l.fq

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	# use syncfilter and syncasm to select organelle reads from HiFi data
	# --stage graph-early \
	if [[ "${_arg_redo}" == "on" ]]; then
		_polap_log3_cmd rm -rf "${_arg_outdir}"
		mkdir -p "${_arg_outdir}"
		touch "${_arg_outdir}/polap.log"
	fi

	local syncasmdir="${_arg_outdir}/syncasm"
	rm -rf "${syncasmdir}"
	mkdir -p "${syncasmdir}"

	# take the input long-read data and save it in outdir for use
	_polap_lib_fastq-normalize-filename "${_arg_long_reads}"
	_arg_long_reads="${_arg_outdir}/tmp/l.fq"
	if [[ ! -s "${_arg_long_reads}" ]]; then
		_polap_log0 "[ASSERT] no file: ${_arg_long_reads}"
		return
	fi

	# annotate the reads
	# get the pt or mt reads
	_polap_lib_readassemble-annotate-read-pt "syncasm"
	_polap_lib_readassemble-get-annotated-read "syncasm"

	local annotatedir="${_arg_outdir}/annotate-read-syncasm"
	seqtk subseq "${_arg_long_reads}" "${annotatedir}/mt.id.all.txt" >"${syncasmdir}/mt.fq"

	# --threads 1 \
	# --threads "${_arg_threads}" \
	local stage="quickview"

	iterate_select_mt() {
		# Iterative recruitment (no separate seed):
		#   select-mt(A_i) -> C_i
		#   gain_i = C_i − A_i
		#   A_{i+1} = A_i ∪ C_i
		#   stop if |gain_i| / |A_{i+1}| < epsilon

		# -------- defaults --------
		local READS="" OUTBASE="" MT_SEED=""
		local EPSILON="0.005" THREADS=8 PRESET="hifi"
		local MAX_ITERS=7
		local PPR_K=121 PPR_S=27 MAX_OCC=200 MIN_SHARED=4 JACC_MIN=0.01
		local NUC_CUT_LOG10=1.30 X_SLOPE=0.20
		local TAIL_NUC=0.05 TAIL_MT=0.30 TAIL_PT=0.30
		local ENSEMBLE_MODE="majority" EMIT_FASTQ=0 HTML_REPORT=""

		local POLAP="${_POLAPLIB_DIR:-.}/polap-bash-select-mt.sh"

		# -------- args --------
		while [[ $# -gt 0 ]]; do
			case "$1" in
			--reads)
				READS="$2"
				shift 2
				;;
			--mt-anchors)
				MT_SEED="$2"
				shift 2
				;; # initial A_0
			--out-prefix)
				OUTBASE="$2"
				shift 2
				;;
			--epsilon)
				EPSILON="$2"
				shift 2
				;;
			--threads)
				THREADS="$2"
				shift 2
				;;
			--preset)
				PRESET="$2"
				shift 2
				;;
			--max-iters)
				MAX_ITERS="$2"
				shift 2
				;;
			--ensemble)
				ENSEMBLE_MODE="$2"
				shift 2
				;;
			--emit-fastq)
				EMIT_FASTQ=1
				shift
				;;
			--html-report)
				HTML_REPORT="$2"
				shift 2
				;;
			--ppr-k)
				PPR_K="$2"
				shift 2
				;;
			--ppr-s)
				PPR_S="$2"
				shift 2
				;;
			--max-occ)
				MAX_OCC="$2"
				shift 2
				;;
			--min-shared)
				MIN_SHARED="$2"
				shift 2
				;;
			--jaccard-min)
				JACC_MIN="$2"
				shift 2
				;;
			--nuc-cut-log10)
				NUC_CUT_LOG10="$2"
				shift 2
				;;
			--x-slope)
				X_SLOPE="$2"
				shift 2
				;;
			--tail-nuc | --x-tail)
				TAIL_NUC="$2"
				shift 2
				;;
			--tail-mt)
				TAIL_MT="$2"
				shift 2
				;;
			--tail-pt)
				TAIL_PT="$2"
				shift 2
				;;
			-h | --help)
				cat <<USAGE
Usage:
  iterate_select_mt --reads R.fq[.gz] --mt seed_ids.txt --out /path/prefix
    [--epsilon 0.005] [--threads 16] [--preset hifi] [--max-iters 30]
    [--ppr-k 121 --ppr-s 27 --max-occ 200 --min-shared 4 --jaccard-min 0.01]
    [--nuc-cut-log10 1.30 --x-slope 0.20]
    [--tail-nuc 0.05] [--tail-mt 0.30] [--tail-pt 0.30]
    [--ensemble majority] [--emit-fastq] [--html-report out/report.html]
USAGE
				return 0
				;;
			*)
				echo "[error] unknown arg: $1" >&2
				return 2
				;;
			esac
		done

		[[ -z "$READS" || -z "$MT_SEED" || -z "$OUTBASE" ]] && {
			echo "[error] require --reads, --mt, --out"
			return 2
		}
		command -v seqtk >/dev/null 2>&1 || {
			echo "[error] need seqtk"
			return 127
		}
		[[ -x "$POLAP" || -f "$POLAP" ]] || {
			echo "[error] not found: $POLAP"
			return 127
		}

		# -------- helpers --------
		uniq_ids() { awk 'NF{print $1}' | tr -d '\r' | LC_ALL=C sort -u; }
		union_ids() { LC_ALL=C sort -u "$@"; }
		diff_ids() { LC_ALL=C comm -23 "$1" "$2"; } # f1 - f2; both sorted
		count_ids() { LC_ALL=C awk 'END{print NR}' "$1"; }

		# -------- init --------
		mkdir -p "$OUTBASE"
		echo "[init] outbase=$OUTBASE  epsilon=$EPSILON  preset=$PRESET  max_iters=$MAX_ITERS"

		local LOGCSV="${OUTBASE}/iter.log.csv"
		if [[ ! -s "$LOGCSV" ]]; then
			echo "iter,a_count,c_count,gain_count,a_next_count,ratio,stopped" >"$LOGCSV"
		fi

		local ITER=0
		local ITERDIR="${OUTBASE}/${ITER}"
		mkdir -p "$ITERDIR"

		# A_0 = uniq(MT_SEED)
		if [[ -s "$MT_SEED" ]]; then
			uniq_ids <"$MT_SEED" >"${ITERDIR}/id.a.txt"
		else
			echo "[error] --mt seed file empty or missing: $MT_SEED" >&2
			return 2
		fi
		: >"${ITERDIR}/id.c.txt"
		: >"${ITERDIR}/id.gain.txt"

		echo "[iter $ITER] |A|=$(count_ids "${ITERDIR}/id.a.txt")"

		# -------- iterate --------
		while [[ $ITER -lt $MAX_ITERS ]]; do
			ITERDIR="${OUTBASE}/${ITER}"
			local NEXTDIR="${OUTBASE}/$((ITER + 1))"
			mkdir -p "$NEXTDIR"

			local A_i="${ITERDIR}/id.a.txt"

			# call select-mt with A_i
			local RUN_PREFIX="${ITERDIR}/run"
			echo "[iter $ITER] select-mt with A_i …"
			bash "${_POLAPLIB_DIR}/polap-bash-select-mt.sh" select-mt \
				--reads "$READS" \
				--mt-anchors "$A_i" \
				--pt-anchors "" \
				--out-prefix "$RUN_PREFIX" \
				--preset "$PRESET" \
				--threads "$THREADS" \
				--ppr-k "$PPR_K" --ppr-s "$PPR_S" \
				--max-occ "$MAX_OCC" --min-shared "$MIN_SHARED" --jaccard-min "$JACC_MIN" \
				--nuc-cut-log10 "$NUC_CUT_LOG10" --x-slope "$X_SLOPE" \
				--x-tail "$TAIL_NUC" \
				--tail-mt "$TAIL_MT" \
				--tail-pt "$TAIL_PT" \
				--ensemble "$ENSEMBLE_MODE" \
				$([[ $EMIT_FASTQ -eq 1 ]] && echo --emit-fastq) \
				${HTML_REPORT:+--html-report "$HTML_REPORT"}

			local C_i="${RUN_PREFIX}.ensemble.mt.ids"
			if [[ ! -s "$C_i" ]]; then
				echo "[warn] no recruits at iter $ITER (C_i empty). stopping."
				# log final
				local aN cN gN aNextN ratio="0.0"
				aN=$(count_ids "$A_i")
				cN=0
				gN=0
				aNextN=$aN
				echo "${ITER},${aN},${cN},${gN},${aNextN},${ratio},1" >>"$LOGCSV"
				break
			fi
			uniq_ids <"$C_i" >"${ITERDIR}/id.c.txt"

			# gain_i = C_i − A_i
			diff_ids "${ITERDIR}/id.c.txt" "${ITERDIR}/id.a.txt" >"${ITERDIR}/id.gain.txt"

			# A_{i+1} = A_i ∪ C_i
			union_ids "${ITERDIR}/id.a.txt" "${ITERDIR}/id.c.txt" >"${NEXTDIR}/id.a.txt"

			# stats & stop
			local aN cN gN aNextN ratio stopped
			aN=$(count_ids "${ITERDIR}/id.a.txt")
			cN=$(count_ids "${ITERDIR}/id.c.txt")
			gN=$(count_ids "${ITERDIR}/id.gain.txt")
			aNextN=$(count_ids "${NEXTDIR}/id.a.txt")

			ratio="0.0"
			if [[ "$aNextN" -gt 0 ]]; then
				ratio=$(echo "scale=6; ${gN} / ${aNextN}" | bc -l)
			fi

			if [[ $(echo "$ratio < $EPSILON" | bc -l) -eq 1 ]]; then
				stopped=1
				echo "[stop] iter=$ITER ratio=$ratio < epsilon=$EPSILON"
			else
				stopped=0
			fi

			echo "[iter $ITER] |C|=$cN  |gain|=$gN  |A_next|=$aNextN  ratio=$ratio"
			echo "${ITER},${aN},${cN},${gN},${aNextN},${ratio},${stopped}" >>"$LOGCSV"

			[[ $stopped -eq 1 ]] && break
			ITER=$((ITER + 1))
		done

		echo "[done] last iter i=$ITER"
		echo "      final A      = ${OUTBASE}/${ITER}/id.a.txt"
		echo "      last C       = ${OUTBASE}/${ITER}/id.c.txt"
		echo "      last gain    = ${OUTBASE}/${ITER}/id.gain.txt"
		echo "      log          = ${OUTBASE}/iter.log.csv"
		seqtk subseq "${_arg_long_reads}" "${OUTBASE}/${ITER}/id.a.txt" | gzip -c >"${OUTBASE}.mt.fastq.gz"
	}

	# Use this if you want to assemble ptDNA using syncfilter method.
	# --pt "${annotatedir}/pt.id.all.txt" \
	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		iterate_select_mt \
			--preset hifi \
			--reads "${_arg_long_reads}" \
			--mt-anchors "${annotatedir}/mt.id.all.txt" \
			--out-prefix "${syncasmdir}/xyz" \
			--threads "${_arg_threads}" \
			--ppr-k 121 --ppr-s 27 --max-occ 200 --min-shared 4 --jaccard-min 0.01 \
			--nuc-cut-log10 1.30 --x-slope 0.20 \
			--x-tail 0.05 \
			--tail-mt 0.30 \
			--tail-pt 0.30 \
			--ensemble majority \
			--emit-fastq \
			--html-report "${syncasmdir}/report.html"
	elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
		# Use this if you want to assemble ptDNA using syncfilter method.
		# --pt "${annotatedir}/pt.id.all.txt" \
		iterate_select_mt \
			--preset ont \
			--reads "${_arg_long_reads}" \
			--mt-anchors "${annotatedir}/mt.id.all.txt" \
			--out-prefix "${syncasmdir}/xyz" \
			--threads "${_arg_threads}" \
			--ppr-k 41 --ppr-s 17 \
			--max-occ 200 --min-shared 4 --jaccard-min 0.01 \
			--nuc-cut-log10 1.30 --x-slope 0.20 \
			--x-tail 0.05 \
			--tail-mt 0.20 \
			--tail-pt 0.20 \
			--ensemble majority \
			--emit-fastq \
			--html-report "${syncasmdir}/report.html"
	else
		_polap_log0 "No such data type: ${_arg_data_type}"
		return
	fi

	# flye on the organelle genome assembly
	# map reads on the seed contigs
	# flye2 on the filtered reads
	rm -rf "${syncasmdir}/mt0"
	# for testing Brassica_napus 1 only
	# _arg_flye_data_type="--nano-hq"
	#
	# Use this if you want to assemble ptDNA using syncfilter method.
	# for organelle in pt mt; do
	for organelle in mt; do
		_polap_log3_cmdout flye "${_arg_flye_data_type}" \
			"${syncasmdir}/xyz.$organelle.fastq.gz" \
			--threads "${_arg_threads}" \
			--out-dir "${syncasmdir}/${organelle}0" \
			--stop-after contigger

		ln -sf syncasm/sample.$stage.$organelle.fastq.gz "${_arg_outdir}/sync.$organelle.fastq.gz"
	done

	################################################################################
	# BEGIN: iteration of assembly
	#
	# input
	local annotatedir="${syncasmdir}"
	# local resolved_fastq="${syncasmdir}/mt.id.all.fq"
	local resolved_fastq="${_arg_long_reads}"
	################################################################################

	local i
	for ((i = 0; i < 6; i++)); do
		local j=$((i + 1))

		# NOTE: annotate for seeding
		# select connected components of the mt contigs only
		_polap_lib_annotate \
			-o "${annotatedir}" \
			-i mt$i

		# NOTE: mito seed
		_polap_lib_seed-mito \
			-o "${annotatedir}" \
			-i mt$i -j mt$j

		if [[ ! -s "${annotatedir}/mt$i/mt.contig.name-mt$j" ]]; then
			_polap_log0 "No mt seed for mt$j"
			break
		fi

		# polap command: assemble-rate
		if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
			_polap_log1 "use then input long reads filtered to remove reads from ptDNA for mtDNA assembly"
			_polap_lib_assemble-rate \
				-o "${annotatedir}" \
				-l "${resolved_fastq}" \
				-w "${_arg_single_min}" \
				-t mt \
				-i mt$i -j mt$j
		elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
			_polap_log1 "use the only selected long reads using organelle gene annotation for mtDNA assembly"
			_polap_lib_assemble-rate \
				-o "${annotatedir}" \
				-l "${annotatedir}"/mt.fq \
				-w "${_arg_single_min}" \
				-t mt \
				-i mt$i -j mt$j
		fi

		if [[ -s "${annotatedir}/mt$j/assembly_graph.gfa" ]]; then
			_polap_lib_mt-extract-dna \
				"${annotatedir}/mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt$j/mtdna"

			_polap_lib_bandage \
				"${annotatedir}/mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt$j/assembly_graph.png"

			ln -sf "mt$j/mtdna/mt.0.fa" \
				"${annotatedir}/mt.$j.fa"

			ln -sf "mt$j/assembly_graph.gfa" \
				"${annotatedir}/mt.$j.gfa"

			ln -sf "mt$j/assembly_graph.png" \
				"${annotatedir}/mt.$j.png"

			_polap_log0 "mtDNA assembly: ${annotatedir}/mt.$j.gfa"
		else
			_polap_log0 "No mt assembly $j"
		fi

	done
	#
	# END: function readassemble-ont-pt-iterate_genus_species
	#######################################################################

	# use the generated reference to assemble a mitochondrial genome.
	local j=$((i + 1))
	# polap command: annotate
	_polap_lib_annotate \
		-o "${annotatedir}" \
		-i mt$i

	# polap command: seed-mito
	_polap_lib_seed-mito \
		-o "${annotatedir}" \
		-i mt$i -j mt$j

	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		# --reference "${annotatedir_pt}"/pt.3.gfa \
		_polap_log1 "use then input long reads filtered to remove reads from ptDNA for mtDNA assembly"
		_polap_lib_assemble-rate \
			-o "${annotatedir}" \
			-l "${resolved_fastq}" \
			-w "${_arg_single_min}" \
			-i mt$i -j mt$j
	elif [[ "${_arg_data_type}" == "nano-raw" ]]; then
		# -l "${_arg_long_reads}" \
		# -l "${annotatedir}"/mt.fq \
		_polap_log1 "use then input long reads with an adjusted omega for the final stage mtDNA assembly"
		_polap_lib_assemble-omega \
			-o "${annotatedir}" \
			-l "${_arg_long_reads}" \
			-i mt$i -j mt$j
	fi

	_polap_lib_annotate \
		-o "${annotatedir}" \
		-i mt$j

	_polap_log0_column "${annotatedir}/mt$j/contig-annotation-depth-table.txt"

	if [[ -s "${annotatedir}/mt$j/assembly_graph.gfa" ]]; then
		_polap_lib_bandage \
			"${annotatedir}/mt$j/assembly_graph.gfa" \
			"${annotatedir}/mt$j/assembly_graph.png"

		ln -sf "mt$j/assembly_graph.gfa" \
			"${annotatedir}/mt.$j.gfa"

		ln -sf "mt$j/assembly_graph.png" \
			"${annotatedir}/mt.$j.png"

		_polap_log0 "mtDNA assembly: ${annotatedir}/mt.$j.gfa"
	else
		_polap_log0 "No MT assembly $j"
	fi

	# extract mtDNA sequence if you can
	#

	# final link
	i=$((j - 1))
	ln -sf "syncasm/mt.$i.gfa" \
		"${_arg_outdir}/mt-sync.0.gfa"

	ln -sf "syncasm/mt.$i.png" \
		"${_arg_outdir}/mt-sync.0.png"

	ln -sf "syncasm/mt.$j.gfa" \
		"${_arg_outdir}/mt-sync.1.gfa"

	ln -sf "syncasm/mt.$j.png" \
		"${_arg_outdir}/mt-sync.1.png"

	# END: iteration of assembly
	################################################################################

	# output
	ln -sf syncasm/sample.$stage.png "${_arg_outdir}/sync.png"
	ln -sf syncasm/pt1/assembly_graph.gfa "${_arg_outdir}/pt-sync.0.gfa"
	ln -sf syncasm/pt1/assembly_graph.fa "${_arg_outdir}/pt-sync.0.fa"
	ln -sf syncasm/mt1/assembly_graph.gfa "${_arg_outdir}/mt-sync.0.gfa"

	_polap_log3_cmd rm -f "${_arg_outdir}/tmp/l.fq"

	bash "${_POLAPLIB_DIR}/polap-bash-select-mt.sh" select-mt-report \
		--out "${syncasmdir}/xyz" \
		--bandage-size 3000x3000 \
		--html-report "${syncasmdir}/xyz.report.html"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
