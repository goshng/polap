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
# Functions for organelle-genome assembly
# Polap executes Flye twice; one for the whole-genome assembly and another for
# the organelle-genome assembly. This script has subcommands for
# organelle-genome assembly. It also has other subcommands that were used
# in the earlier versions.
#
# function _polap_oga_determine-long-read-file {
# function _run_polap_prepare-seeds { # prepare seed contigs in a not usual way
# function _run_polap_map-reads { # selects reads mapped on a genome assembly
# function _run_polap_test-reads { # selects reads mapped on a genome assembly
# function _run_polap_select-reads { # selects reads mapped on a genome assembly
# function _run_polap_flye2 { # executes Flye for an organelle-genome assembly
# function _run_polap_flye-polishing { # finish a Flye organelle-genome assembly upto flye-polishing step
# function _run_polap_report-assembly { # report an organelle-genome assembly result
# function _run_polap_x-v0.3.7-collect-reads() { # replaced by select-reads
# function _run_polap_x-v0.2.6-select-reads() { # selects reads mapped on a genome assembly
# function _run_polap_x-select-reads { # selects reads mapped on a genome assembly
# function _run_polap_x-v0.2.6-flye2() { # executes Flye for an organelle-genome assembly
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

source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"

function _polap_oga_determine-long-read-file {
	local -n result_ref=$1

	if [[ "${_arg_long_reads_is}" == "off" ]]; then
		if [[ -s "${_polap_var_outdir_lk_fq_gz}" ]]; then
			_polap_log2 "    we utilize all available size-limited long-read data for our analysis: ${_polap_var_outdir_lk_fq_gz}"
			result_ref="${_polap_var_outdir_lk_fq_gz}"
		elif [[ -s "${_polap_var_outdir_nk_fq_gz}" ]]; then
			_polap_log2 "    we utilize the sampled and size-limited long-read data for our analysis: ${_polap_var_outdir_nk_fq_gz}"
			result_ref="${_polap_var_outdir_nk_fq_gz}"
		else
			die "ERROR: no such file: ${_polap_var_outdir_lk_fq_gz}, ${_polap_var_outdir_nk_fq_gz}"
		fi
	else
		_polap_log2 "    we utilize the long-read data supplied through the command-line option -l."
		result_ref="${_arg_long_reads}"
	fi
}

################################################################################
# Tests a few read-selection methods.
#
# 1. outdir/1/01-reads/1/3000/ upto a large value (3, 5, 10 values)
#      ptGAUL with different V11 starting from 3k
#      V10/V11 > 0.7 might need to be tested as well: 0.5~0.9 for 5 values
# 2. outdir/1/01-reads/2/3000/
#      POLAP intra-contig with different V11 starting from 3k
#      V10/V11 > 0.7 might need to be tested as well
#      (V4-V3)/V2 > 0.7 might need to be tested as well
# 3. outdir/1/01-reads/3/3000/
#      POLAP inter-contig with different V11 starting from 3k
#      V10/V11 > 0.7 might need to be tested as well
# 4. outdir/1/01-reads/4/3000/
#      POLAP inter-contig with different V2 starting from 3k
#      V10/V11 > 0.7 might need to be tested as well
#
# View:
# reads stats: number, range of lengths, mean, median
# assembly: total bases, number of fragments
#
# Arguments:
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   --
#   -r ${_arg_pair_min}: minimum minimap2 alignment length for a pair of contigs
#   -x ${_arg_bridge_min}: minimum long-read length for connecting the pair of contigs
#   -w ${_arg_single_min}: minimum minimap2 alignment length for a single contig
#
# Inputs:
#   _polap_var_oga_reads="${_polap_var_oga}/02-reads"
#   ${_polap_var_oga_reads}/contig.tab
#
# Outputs:
#   ${MTSEEDSDIR}
#   ${MTDIR}/contig.fa
#   ${MTDIR}/contig_total_length.txt
#   ${MTDIR}/contig.paf
#   ${MTDIR}/contig.tab
#   ${MTSEEDSDIR}/1.names
#   ${MTSEEDSDIR}/1.fq.gz
#   ${MTSEEDSDIR}/2.fq.gz
#
# Outputs:
#   ${_polap_var_oga_reads}/contig.tab
################################################################################
function _run_polap_x-v0.3.7-collect-reads() { # replaced by select-reads
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	# for contigs
	#	_polap_var_ga_contigger_edges_fasta=o/30-contigger/contigs.fasta
	# for edges

	help_message=$(
		cat <<HEREDOC
# Select long reads mapped on a genome assembly.
#
# Arguments:
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number (or 0)
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   -l ${_arg_long_reads}
# Inputs:
#   ${_polap_var_mtcontigname}
#   ${_polap_var_ga_contigger_edges_fasta}
#   input long read data: 1. ${_polap_var_outdir_lk_fq_gz} <- input data used for the whole-genome assembly
#                         2. ${_arg_long_reads}          <- input long-read data
# Outputs:
#   ${_polap_var_oga_seeds}/ptgaul.fq.gz   <- ptGAUL read-selection
#   ${_polap_var_oga_seeds}/single.fq.gz   <- single-mapped reads, if we have many such reads
#   ${_polap_var_oga_seeds}/pair.fq.gz     <- only pair-mapped reads (not useful)
#   ${_polap_var_oga_seeds}/combined.fq.gz <- single-mapped reads plus pair-mapped reads
Example: $(basename "$0") ${_arg_menu[0]} -i ${_arg_inum} -j ${_arg_jnum}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	if [[ "${_arg_menu[1]}" == "view" ]]; then
		_polap_log0 "view: ptGAUL mapped reads"
		if [[ -s "${_polap_var_oga_seeds}/ptgaul.fq.gz" ]]; then
			seqkit stats "${_polap_var_oga_seeds}/ptgaul.fq.gz"
		else
			_polap_log0 "No such file: ${_polap_var_oga_seeds}/ptgaul.fq.gz"
		fi
		_polap_log0 "view: single-mapped reads"
		if [[ -s "${_polap_var_oga_seeds}/single.fq.gz" ]]; then
			seqkit stats "${_polap_var_oga_seeds}/single.fq.gz"
		else
			_polap_log0 "No such file: ${_polap_var_oga_seeds}/single.fq.gz"
		fi
		_polap_log0 "view: pair-mapped reads"
		if [[ -s "${_polap_var_oga_seeds}/pair.fq.gz" ]]; then
			seqkit stats "${_polap_var_oga_seeds}/pair.fq.gz"
		else
			_polap_log0 "No such file: ${_polap_var_oga_seeds}/pair.fq.gz"
		fi
		_polap_log0 "view: combined-mapped reads"
		if [[ -s "${_polap_var_oga_seeds}/combined.fq.gz" ]]; then
			seqkit stats "${_polap_var_oga_seeds}/combined.fq.gz"
		else
			_polap_log0 "No such file: ${_polap_var_oga_seeds}/combined.fq.gz"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	_polap_log0 "selecting long-reads mapped on the seed contigs in file"

	_polap_log1 "  please, wait for a long-read data selection ..."
	_polap_log1 "    ${_arg_inum} (source) -> ${_arg_jnum} (target) ..."
	_polap_log1 "    bridge=${_arg_bridge_min}"
	_polap_log1 "    p_mapping=${_arg_pair_min}"
	_polap_log1 "    s_mapping=${_arg_single_min}"
	_polap_log1 "    min_len_read=${_arg_min_read_length}"

	# MT: MPAIR=3000 MBRIDGE=3000 MSINGLE=3000
	# PT: MPAIR=1000 MBRIDGE=5000 MSINGLE=0
	_polap_log1 "  selecting seeds that meet the selection criteria ..."
	_polap_log2 "    bridge (POLAP pairs bridge minimum)=${_arg_bridge_min}"
	_polap_log2 "    p_mapping (POLAP pairs alignment minimum)=${_arg_pair_min}"
	_polap_log2 "    s_mapping (POLAP single alignment minimum)=${_arg_single_min}"
	_polap_log2 "    min_len_read=${_arg_min_read_length}"
	_polap_log2 "    input1: ${_polap_var_mtcontigname}"
	_polap_log2 "    input2: ${_polap_var_oga_reads}/contig.tab"
	_polap_log2 "    output1: ${_polap_var_oga_reads}/ptgaul.names"
	_polap_log2 "    output2: ${_polap_var_oga_reads}/single.names for reads mapped on a single contig"
	_polap_log2 "    output3: ${_polap_var_oga_reads}/pair.names"
	_polap_log2 "    output4: ${_polap_var_oga_reads}/combined.names"
	_polap_log2 "    output5: ${_polap_var_oga_reads}/<edge1>-<edge2>.name for reads mapped on contig pairs"
	_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-pairs.R \
		-m ${_polap_var_mtcontigname} \
		-t ${_polap_var_oga_reads}/contig.tab \
		-o ${_polap_var_oga_reads} \
		-r ${_arg_pair_min} \
		-x ${_arg_bridge_min} \
    -w ${_arg_single_min} \
    --all \
    >${_polap_output_dest} 2>&1"

	_polap_log0 "extracting reads using read names ..."

	if [[ "${_arg_long_reads_is}" == "off" ]]; then
		local _source_long_reads_fq="${_polap_var_outdir_lk_fq_gz}"
	else
		local _source_long_reads_fq="${_arg_long_reads}"
	fi

	if [[ -s "${_source_long_reads_fq}" ]]; then
		_polap_log1 "  using ${_source_long_reads_fq} for a source of long-read data"
	else
		_polap_log0 "ERROR: no such file: ${_source_long_reads_fq}"
		return
	fi

	_polap_log3_cmd mkdir -p "${_polap_var_oga_seeds}"

	_polap_log1 "  extracting the single-mapped ..."
	_polap_log2 "    input1: ${_source_long_reads_fq}"
	_polap_log2 "    input2: ${_polap_var_oga_reads}/ptgaul.names"
	_polap_log2 "    output: ${_polap_var_oga_seeds}/ptgaul.fq.gz"

	_polap_log3_pipe "seqtk subseq \
    ${_source_long_reads_fq} \
    ${_polap_var_oga_reads}/ptgaul.names |\
    gzip >${_polap_var_oga_seeds}/ptgaul.fq.gz"

	_polap_log1 "  extracting the single-mapped ..."
	_polap_log2 "    input1: ${_source_long_reads_fq}"
	_polap_log2 "    input2: ${_polap_var_oga_reads}/single.names"
	_polap_log2 "    output: ${_polap_var_oga_seeds}/single.fq.gz"

	_polap_log3_pipe "seqtk subseq \
    ${_source_long_reads_fq} \
    ${_polap_var_oga_reads}/single.names |\
    gzip >${_polap_var_oga_seeds}/single.fq.gz"

	_polap_log1 "  extracting the single-mapped and pair-mapped reads ..."
	_polap_log2 "    input1: ${_source_long_reads_fq}"
	_polap_log2 "    input2: ${_polap_var_oga_reads}/combined.names"
	_polap_log2 "    output: ${_polap_var_oga_seeds}/combined.fq.gz"

	_polap_log3_pipe "seqtk subseq \
    ${_source_long_reads_fq} \
    ${_polap_var_oga_reads}/combined.names |\
    gzip >${_polap_var_oga_seeds}/combined.fq.gz"

	_polap_log1 NEXT: $(basename "$0") flye2 -o "${_arg_outdir}" -j "${_arg_jnum}" -t "${_arg_threads}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Selects reads mapped on a genome assembly.
#
# Arguments:
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   -r ${_arg_pair_min}: minimum minimap2 alignment length for a pair of contigs
#   -x ${_arg_bridge_min}: minimum long-read length for connecting the pair of contigs
#   -w ${_arg_single_min}: minimum minimap2 alignment length for a single contig
# Inputs:
#   ${_polap_var_mtcontigname}
#   ${_polap_var_ga_contigger_edges_fasta}
# Outputs:
#   ${MTSEEDSDIR}
#   ${MTDIR}/contig.fa
#   ${MTDIR}/contig_total_length.txt
#   ${MTDIR}/contig.paf
#   ${MTDIR}/contig.tab
#   ${MTSEEDSDIR}/1.names
#   ${MTSEEDSDIR}/1.fq.gz
#   ${MTSEEDSDIR}/2.fq.gz
################################################################################
function _run_polap_x-v0.2.6-select-reads() { # selects reads mapped on a genome assembly
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	local MTDIR="${_polap_var_oga}"                                              # target: ${_polap_var_oga}
	local MTSEEDSDIR="${_polap_var_oga}/seeds"                                   # ${_polap_var_seeds} for oga-class
	local _polap_var_mtcontigname="${_polap_var_ga}/mt.contig.name-${_arg_jnum}" # ${_polap_var_mtcontigname}

	# for contigs
	#	_polap_var_ga_contigger_edges_fasta=o/30-contigger/contigs.fasta
	# for edges

	help_message=$(
		cat <<HEREDOC
# Select long reads mapped on a genome assembly.
#
# Arguments:
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number (or 0)
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   -w ${_arg_single_min}: set the option values of -r, -x, -w to the same one
# Inputs:
#   ${_polap_var_mtcontigname}
#   ${_polap_var_ga_contigger_edges_fasta}
# Outputs:
#   ${MTDIR}/contig.fa
#   ${MTSEEDSDIR}/1.fq.gz <- single-mapped reads
#   ${MTSEEDSDIR}/2.fq.gz <- reduced single-mapped reads, if we have many such reads
#   ${MTSEEDSDIR}/3.fq.gz <- reduced single-mapped reads plus pair-mapped reads
Example: $(basename "$0") ${_arg_menu[0]} -i ${_arg_inum} -j ${_arg_jnum} -w ${_arg_single_min}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	_polap_log0 "selecting long-reads mapped on the seed contigs in file: ${_polap_var_mtcontigname} and ${_polap_var_ga_contigger_edges_fasta}"
	_polap_log1 "  input1: ${_polap_var_mtcontigname}"
	_polap_log1 "  input2: ${_polap_var_ga_contigger_edges_fasta}"

	if [ ! -s "${_polap_var_mtcontigname}" ]; then
		_polap_log0 "ERROR: no such mt.contig.name file: ${_polap_var_mtcontigname}"
		exit $EXIT_ERROR
	fi

	if [ ! -s "${_polap_var_ga_contigger_edges_fasta}" ]; then
		_polap_log0 "ERROR: no assembly fasta file: ${_polap_var_ga_contigger_edges_fasta}"
		exit $EXIT_ERROR
	fi

	# if [ -d "${MTDIR}" ] && [ "${_arg_yes}" = "off" ]; then
	if [ -d "${MTDIR}" ] && [ "${_arg_redo}" = "off" ]; then
		while true; do
			read -r -p "Folder [${MTDIR}] already exists. Do you want to replace it? [y/n] " yn
			case $yn in
			[Yy]*)
				_polap_log3_cmd rm -rf "${MTDIR}"
				_polap_log0 "  folder ${MTDIR} is deleted."
				break
				;;
			[Nn]*)
				_polap_log0 "  folder ${MTDIR} is not replaced."
				_polap_log0 "  you might want a new mt.contig.name file for flye2 step."
				exit $EXIT_FAIL
				;;
			*) _polap_log0 "Please answer yes or no." ;;
			esac
		done
	else
		_polap_log3_cmd rm -rf "${MTDIR}"
		_polap_log1 "  ${MTDIR} is deleted if there is one."
	fi

	_polap_log1 "  creating ${_arg_outdir}/organelle-assembly_${_arg_inum}-${_arg_jnum}"
	echo "$CMD" >"${_arg_outdir}/organelle-assembly_${_arg_inum}-${_arg_jnum}"

	_polap_log1 "  creates ${MTSEEDSDIR}"
	_polap_log3_cmd mkdir -p "${MTSEEDSDIR}"
	_polap_log3_cmd ln -s "$PWD/${_polap_var_outdir_nk_fq_gz}" -t "${MTDIR}"
	_polap_log1 "  extracts contig sequeces from the assembly: ${_polap_var_ga_contigger_edges_fasta}"
	_polap_log2 "    input1: ${_polap_var_ga_contigger_edges_fasta}"
	_polap_log2 "    input2: ${_polap_var_mtcontigname}"
	_polap_log2 "    output: ${MTDIR}/contig.fa"
	_polap_log3_pipe "seqkit grep \
    --threads ${_arg_threads} \
    -f ${_polap_var_mtcontigname} \
		${_polap_var_ga_contigger_edges_fasta} \
		-o ${MTDIR}/contig.fa \
		2>${_polap_output_dest}"

	# we could circularize a single contig.
	local contig_count=$(wc -l <"${_polap_var_mtcontigname}")
	_polap_log1 "  number of seed contigs: ${contig_count}"
	if [[ ${_arg_circularize} == "on" ]]; then
		if [ "$contig_count" -eq 1 ]; then

			_polap_log1 "  circularizes the single seed contig ..."
			_polap_log2 "    input1: ${MTDIR}"
			_polap_log2 "    input2: ${MTDIR}/contig.fa"
			_polap_log2 "    input3: ${_polap_var_mtcontigname}"
			_polap_log2 "    output1: ${_polap_var_mtcontigname}-backup"
			_polap_log2 "    output2: (new) ${_polap_var_mtcontigname}"
			_polap_log2 "    output3: ${MTDIR}/contig.fa-backup"
			_polap_log2 "    output4: (new) ${MTDIR}/contig.fa"
			_polap_log3_cmd bash "${_POLAPLIB_DIR}/polap-bash-half-cut.sh" "${MTDIR}" "${MTDIR}/contig.fa" "${_polap_var_mtcontigname}"

			# cp "${MTDIR}/contig.fa" "${MTDIR}/contig.fa-backup"
			# seqkit fx2tab --length --name "${MTDIR}"/contig.fa -o "${MTDIR}"/contig.fa.len >/dev/null 2>&1
			# A=$(cut -f2 "${MTDIR}"/contig.fa.len)
			# B=$(echo "scale=0; $A/2" | bc)
			# C=$((B + 1))
			# seqkit subseq -r 1:"$B" "${MTDIR}"/contig.fa -o "${MTDIR}"/c1.fa >/dev/null 2>&1
			# seqkit subseq -r "$C":"$A" "${MTDIR}"/contig.fa -o "${MTDIR}"/c2.fa >/dev/null 2>&1
			# cat "${MTDIR}"/c?.fa | seqkit replace -p '.+' -r 'edge_{nr}' -o "${MTDIR}"/contig.fa >/dev/null 2>&1
			# cp "${_polap_var_mtcontigname}" "${_polap_var_mtcontigname}"-backup
			# echo -e "edge_1\nedge_2" >"${_polap_var_mtcontigname}"

			_polap_log1 "    creating new ${MTDIR}/contig.fa and ${_polap_var_mtcontigname}"

		else
			_polap_log0 "Not implemented yet!"
			exit $EXIT_ERROR
			# "${_POLAPLIB_DIR}"/run-polap-single.R "${MTSEEDSDIR}"/contig.tab "${MTSEEDSDIR}" "${_arg_single_min}" >/dev/null 2>&1
			# cat "${MTSEEDSDIR}"/single.names | sort | uniq >"${MTSEEDSDIR}"/1.names
			# echo "INFO: creates long read single name in ${MTSEEDSDIR}/1.names"
		fi
	fi

	_polap_log1 "  please, wait for a long-read data selection ..."
	_polap_log1 "    ${_arg_inum} (source) -> ${_arg_jnum} (target) ..."
	_polap_log1 "    bridge=${_arg_bridge_min}"
	_polap_log1 "    p_mapping=${_arg_pair_min}"
	_polap_log1 "    s_mapping=${_arg_single_min}"
	_polap_log1 "    min_len_read=${_arg_min_read_length}"

	_polap_utility_get_contig_length \
		"${MTDIR}/contig.fa" \
		"${MTDIR}/contig_total_length.txt"
	local CONTIG_LENGTH=$(<"${MTDIR}/contig_total_length.txt")
	_polap_log2_cat "${_polap_var_oga}/contig_total_length.txt"
	# CONTIG_LENGTH=$(_polap_utility_get_contig_length "${MTDIR}/contig.fa")
	# echo "$CONTIG_LENGTH" >"${MTDIR}"/contig_total_length.txt
	_polap_log1 "  organelle genome size based on the seed contig selection: $CONTIG_LENGTH"

	_polap_log1 "  mapping long-read data on the seed contigs using minimap2 ..."
	_polap_log2 "    input1: ${MTDIR}/contig.fa"
	_polap_log2 "    input2: ${_polap_var_outdir_nk_fq_gz}"
	_polap_log2 "    output: ${MTDIR}/contig.paf"
	if [[ -s "${MTDIR}"/contig.paf ]] && [[ "${_arg_redo}" = "off" ]]; then
		_polap_log1 "  found: ${MTDIR}/contig.paf, skipping the minimap2 mapping step ..."
	else
		_polap_log3_pipe "minimap2 -cx \
      ${_arg_minimap2_data_type} \
      ${MTDIR}/contig.fa \
      ${_polap_var_outdir_nk_fq_gz} \
      -t ${_arg_threads} \
      -o ${MTDIR}/contig.paf \
      >${_polap_output_dest} 2>&1"
	fi

	_polap_log1 "  filtering out reads less than ${_arg_min_read_length} bp and converting PAF to TAB ..."
	_polap_log2 "    input1: ${MTDIR}/contig.paf"
	_polap_log2 "    output: ${MTDIR}/contig.tab"
	# cut -f1-11 "${MTDIR}"/contig.paf | awk -v minlength="${_arg_min_read_length}" '{if ($2>=minlength) {print}}' >"${MTDIR}"/contig.tab
	_polap_log3_cmd bash "${_POLAPLIB_DIR}/polap-bash-minimap2-paf2tab.sh" "${_arg_min_read_length}" "${MTDIR}/contig.paf" "${MTDIR}/contig.tab"

	# MT: MPAIR=3000 MBRIDGE=3000 MSINGLE=3000
	# PT: MPAIR=1000 MBRIDGE=5000 MSINGLE=0
	_polap_log1 "  selecting seeds that meet the selection criteria ..."
	_polap_log2 "    bridge (POLAP pairs bridge minimum)=${_arg_bridge_min}"
	_polap_log2 "    p_mapping (POLAP pairs alignment minimum)=${_arg_pair_min}"
	_polap_log2 "    s_mapping (POLAP single alignment minimum)=${_arg_single_min}"
	_polap_log2 "    min_len_read=${_arg_min_read_length}"
	_polap_log2 "    input1: ${_polap_var_mtcontigname}"
	_polap_log2 "    input2: ${MTDIR}/contig.tab"
	_polap_log2 "    output: ${MTSEEDSDIR}/single.names for reads mapped on a single contig"
	_polap_log2 "    output: ${MTSEEDSDIR}/<edge1>-<edge2>.name for reads mapped on contig pairs"
	_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-pairs.R \
		${_polap_var_mtcontigname} \
		${MTDIR}/contig.tab \
		${MTSEEDSDIR} \
		${_arg_pair_min} \
    ${_arg_bridge_min} \
    ${_arg_single_min} \
    >${_polap_output_dest} 2>&1"
	# "${_POLAPLIB_DIR}"/run-polap-pairs.R "${_polap_var_mtcontigname}" ${MTDIR}/contig.tab ${MTSEEDSDIR} ${_arg_pair_min} ${_arg_bridge_min} ${_arg_single_min} >/dev/null 2>&1

	# cat "${MTSEEDSDIR}/"*".name" "${MTSEEDSDIR}"/single.names | sort | uniq >"${MTSEEDSDIR}"/1.names
	_polap_log1 "  creates names of reads mapped on a single name in ${MTSEEDSDIR}/1.names"
	_polap_log2 "    input1: ${MTSEEDSDIR}/single.names"
	_polap_log2 "    output: ${MTSEEDSDIR}/1.names"
	_polap_log3_pipe "cat ${MTSEEDSDIR}/single.names | sort | uniq >${MTSEEDSDIR}/1.names"

	# seqkit grep --threads ${_arg_threads} -f "${MTSEEDSDIR}"/1.names ${_polap_var_outdir_nk_fq_gz} -o "${MTSEEDSDIR}"/1.fq.gz >/dev/null 2>&1
	_polap_log1 "  creating reads mapped on a single contig in ${MTSEEDSDIR}/1.fq.gz"
	_polap_log2 "    input1: ${_polap_var_outdir_nk_fq_gz}"
	_polap_log2 "    input2: ${MTSEEDSDIR}/1.names"
	_polap_log2 "    output: ${MTSEEDSDIR}/1.fq.gz"
	seqtk subseq "${_polap_var_outdir_nk_fq_gz}" "${MTSEEDSDIR}"/1.names | gzip >"${MTSEEDSDIR}"/1.fq.gz

	_polap_log1 "  computing the total size of the single-mapped reads ..."
	_polap_log2 "    input1: ${MTSEEDSDIR}/1.fq.gz"
	local _contig_length_bp=$(_polap_utility_convert_bp ${CONTIG_LENGTH})
	_polap_log2 "    input2 (seed contig size): ${_contig_length_bp}"
	local TOTAL_LENGTH=$(seqkit stats -Ta "${MTSEEDSDIR}"/1.fq.gz | csvtk cut -t -f "sum_len" | csvtk del-header)
	local _total_length_bp=$(_polap_utility_convert_bp ${TOTAL_LENGTH})
	_polap_log2 "    result1 (total size of reads mapped on a single contig): ${_total_length_bp}"
	local EXPECTED_ORGANELLE_COVERAGE=$((TOTAL_LENGTH / CONTIG_LENGTH))
	_polap_log2 "    result2 (expected organelle coverage): ${EXPECTED_ORGANELLE_COVERAGE}x"

	_polap_log1 "  reducing the single-mapped reads upto the coverage of ${_arg_coverage_oga}x"
	_polap_log2 "    input1: ${MTSEEDSDIR}/1.fq.gz"
	_polap_log2 "    output: ${MTSEEDSDIR}/2.fq.gz"
	if [[ "${_arg_test}" == "on" ]]; then
		_polap_log0 "    OPTION: --test : No reduction of the test long-read data"
		_polap_log3_cmd ln -s $(realpath "${MTSEEDSDIR}"/1.fq.gz) "${MTSEEDSDIR}"/2.fq.gz
	elif [[ "${_arg_coverage_check}" == "off" ]]; then
		_polap_log0 "    OPTION: --no-coverage-check : No reduction of the long-read data"
		_polap_log3_cmd ln -s $(realpath "${MTSEEDSDIR}"/1.fq.gz) "${MTSEEDSDIR}"/2.fq.gz
	else
		if [ "$EXPECTED_ORGANELLE_COVERAGE" -lt "${_arg_coverage_oga}" ]; then
			_polap_log1 "    no reduction of the long-read data because $EXPECTED_ORGANELLE_COVERAGE < ${_arg_coverage_oga}"
			_polap_log3_cmd ln -s $(realpath "${MTSEEDSDIR}"/1.fq.gz) "${MTSEEDSDIR}"/2.fq.gz
		else
			_polap_log1 "    SUGGESTION: you might want to increase the minimum read lengths (use -w or -m) because you have enough long-read data."
			RATE=$(echo "scale=10; ${_arg_coverage_oga}/$EXPECTED_ORGANELLE_COVERAGE" | bc)
			_polap_log0 "  long-read data reduction by rate of $RATE <= COV[${_arg_coverage_oga}] / long-read organelle coverage[$EXPECTED_ORGANELLE_COVERAGE]"
			_polap_log1 "    sampling long-read data by $RATE ... wait ..."
			# seqkit sample -p "$RATE" "${MTSEEDSDIR}/1.fq.gz" -o "${MTSEEDSDIR}/2.fq.gz" >/dev/null 2>&1
			_polap_lib_random-get
			local seed=${_polap_var_random_number}
			_polap_log1 "    random seed for reducing the single-mapped long-read data: ${seed}"
			rm -f "${MTSEEDSDIR}/2.fq.gz" 
			_polap_log3_pipe "seqkit sample \
        -p ${RATE} \
        -s ${seed} \
        ${MTSEEDSDIR}/1.fq.gz \
        -o ${MTSEEDSDIR}/2.fq.gz 2>${_polap_output_dest}"
		fi
	fi

	# 1.fq.gz: reads mapped on a single read
	# 2.fq.gz: a reduced data of 1.fq.gz
	#          or a reduced data of 1.fq.gz plus the inter-mapped reads, if any
	# single.names.2: read names from the reduced data
	#
	# v0.2.6
	# input1: single.names: single-mapped read names
	# input2: edge-number1_edge_number2.name: pair-mapped read names
	# 1.names: unique read names from single.names
	# 1.fq.gz: reads using 1.names
	# 2.fq.gz: reduced from 1.fq.gz but later replaced by some with additional pair-mapped reads
	# single.names.2: names from 2.fq.gz with only single-mapped reads or reduced names
	# 1.names.2: single.names.2 + all pair-mapped read names
	#
	# single.names.2: single-mapped reads -> read data
	# 1.names.2: reads used for organelle-genome
	# use:
	# seqtk subseq ${_polap_var_outdir_nk_fq_gz} ${MTSEEDSDIR}/1.names.2 | gzip >${MTSEEDSDIR}/3.fq.gz
	_polap_log1 "  adding pair-mapped bridging long reads to the single-mapped data ..."
	local C=$(ls -1 "${MTSEEDSDIR}/"*".name" 2>/dev/null | wc -l)
	if [ "$C" != 0 ]; then
		_polap_log2 "    input1: combinations of $C bridging reads"
		_polap_log1 "  creating a list of reduced reads mapped on a single contig"
		_polap_log2 "    input1: ${MTSEEDSDIR}/2.fq.gz"
		_polap_log2 "    output: ${MTSEEDSDIR}/single.names.2"
		_polap_log3_cmd seqkit seq -n -i "${MTSEEDSDIR}"/2.fq.gz -o "${MTSEEDSDIR}"/single.names.2

		_polap_log3 "  FIXME: _polap_log3_cmd or _polap_log3_pipe for cat *"
		cat "${MTSEEDSDIR}/"*".name" "${MTSEEDSDIR}"/single.names.2 | sort | uniq >"${MTSEEDSDIR}/1.names.2"
		# seqkit grep --threads ${_arg_threads} -f "${MTSEEDSDIR}"/1.names.2 ${_polap_var_outdir_nk_fq_gz} -o "${MTSEEDSDIR}"/2.fq.gz >/dev/null 2>&1

		_polap_log1 "  extracting single- and pair-mapped reads ..."
		_polap_log2 "    input1: ${_polap_var_outdir_nk_fq_gz}"
		_polap_log2 "    input1: ${MTSEEDSDIR}/1.names.2"
		_polap_log2 "    output: ${MTSEEDSDIR}/3.fq.gz"
		_polap_log3_pipe "seqtk subseq ${_polap_var_outdir_nk_fq_gz} ${MTSEEDSDIR}/1.names.2 | gzip >${MTSEEDSDIR}/3.fq.gz"
	else
		_polap_log1 "    no pair-mapped reads: 2.fq.gz -> 3.fq.gz"
		_polap_log3_cmd ln -s $(realpath "${MTSEEDSDIR}"/2.fq.gz) "${MTSEEDSDIR}"/3.fq.gz
	fi
	_polap_log1 "  organelle reads in ${MTSEEDSDIR}/3.fq.gz"

	# put the backup to the original
	if [[ ${_arg_circularize} == "on" ]]; then
		_polap_log1 "  circularize option: putting seed contigs back to the original mt.contig.name and ${MTDIR}/contig.fa"
		if [[ -s "${_polap_var_mtcontigname}"-backup ]]; then
			mv "${_polap_var_mtcontigname}"-backup "${_polap_var_mtcontigname}"
			mv "${MTDIR}/contig.fa-backup" "${MTDIR}/contig.fa"
		else
			echo "DEV: not implemented yet"
			exit $EXIT_ERROR
		fi
	fi

	_polap_log1 NEXT: $(basename "$0") flye2 -o "${_arg_outdir}" -j "${_arg_jnum}" -t "${_arg_threads}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# simply use read names to extract reads from a long-read data file
################################################################################
function _run_polap_x-select-reads { # selects reads mapped on a genome assembly
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	local MTDIR="${_polap_var_oga}"                                              # target: ${_polap_var_oga}
	local MTSEEDSDIR="${_polap_var_oga}/seeds"                                   # ${_polap_var_seeds} for oga-class
	local _polap_var_mtcontigname="${_polap_var_ga}/mt.contig.name-${_arg_jnum}" # ${_polap_var_mtcontigname}

	# for contigs
	#	_polap_var_ga_contigger_edges_fasta=o/30-contigger/contigs.fasta
	# for edges

	help_message=$(
		cat <<HEREDOC
# Extract long reads mapped on a genome assembly.
#
# Arguments:
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   -l ${_arg_long_reads}
# Inputs:
#   input long read data: 1. ${_polap_var_outdir_nk_fq_gz} <- input data used for the whole-genome assembly
#                         2. ${_arg_long_reads}          <- input long-read data
#   ${_arg_outdir}/${_arg_jnum}/seeds/1.names         <- single-mapped read
#   ${_arg_outdir}/${_arg_jnum}/seeds/single.names.2  <- reduced single-mapped read
#   ${_arg_outdir}/${_arg_jnum}/seeds/1.names.2       <- reduced single-mapped read + pair-mapped read
#   ${_polap_var_ga_contigger_edges_fasta}
# Outputs:
#   ${_polap_var_oga_seeds}/2.fq.gz <- reduced single-mapped reads, if we have many such reads
#   ${_polap_var_oga_seeds}/3.fq.gz <- reduced single-mapped reads plus pair-mapped reads
Example: $(basename "$0") ${_arg_menu[0]} -j ${_arg_jnum}
Example: $(basename "$0") ${_arg_menu[0]} -j ${_arg_jnum} -l l.fq
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	if [[ "${_arg_menu[1]}" == "view" ]]; then
		_polap_log0 "view: single-mapped reads"
		if [[ -s "${_polap_var_oga_seeds}/1.fq.gz" ]]; then
			seqkit stats "${_polap_var_oga_seeds}/1.fq.gz"
		else
			_polap_log0 "No such file: ${_polap_var_oga_seeds}/1.fq.gz"
		fi
		_polap_log0 "view: reduced single-mapped reads"
		if [[ -s "${_polap_var_oga_seeds}/2.fq.gz" ]]; then
			seqkit stats "${_polap_var_oga_seeds}/2.fq.gz"
		else
			_polap_log0 "No such file: ${_polap_var_oga_seeds}/2.fq.gz"
		fi
		_polap_log0 "view: reduced single-mapped + pair-mapped reads"
		if [[ -s "${_polap_var_oga_seeds}/3.fq.gz" ]]; then
			seqkit stats "${_polap_var_oga_seeds}/3.fq.gz"
		else
			_polap_log0 "No such file: ${_polap_var_oga_seeds}/3.fq.gz"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	_polap_log0 "extracting reads using read names ..."

	if [[ "${_arg_long_reads_is}" == "off" ]]; then
		local _source_long_reads_fq="${_polap_var_outdir_nk_fq_gz}"
	else
		local _source_long_reads_fq="${_arg_long_reads}"
	fi
	if [[ -s "${_source_long_reads_fq}" ]]; then
		_polap_log1 "  using ${_source_long_reads_fq} for a source of long-read data"
	else
		_polap_log0 "ERROR: no such file: ${_source_long_reads_fq}"
		return
	fi

	_polap_log3_cmd mkdir -p "${_polap_var_oga_seeds}"

	_polap_log1 "  extracting the single-mapped ..."
	_polap_log2 "    input1: ${_source_long_reads_fq}"
	_polap_log2 "    input2: ${MTSEEDSDIR}/1.names"
	_polap_log2 "    output: ${_polap_var_oga_seeds}/1.fq.gz"

	_polap_log3_pipe "seqtk subseq \
    ${_source_long_reads_fq} \
    ${MTSEEDSDIR}/1.names |\
    gzip >${_polap_var_oga_seeds}/1.fq.gz"

	_polap_log1 "  extracting the reduced single-mapped ..."
	_polap_log2 "    input1: ${_source_long_reads_fq}"
	_polap_log2 "    input2: ${MTSEEDSDIR}/single.names.2"
	_polap_log2 "    output: ${_polap_var_oga_seeds}/2.fq.gz"

	_polap_log3_pipe "seqtk subseq \
    ${_source_long_reads_fq} \
    ${MTSEEDSDIR}/single.names.2 |\
    gzip >${_polap_var_oga_seeds}/2.fq.gz"

	_polap_log1 "  extracting the reduced single-mapped and pair-mapped reads ..."
	_polap_log2 "    input1: ${_source_long_reads_fq}"
	_polap_log2 "    input2: ${MTSEEDSDIR}/1.names.2"
	_polap_log2 "    output: ${_polap_var_oga_seeds}/3.fq.gz"

	_polap_log3_pipe "seqtk subseq \
    ${_source_long_reads_fq} \
    ${MTSEEDSDIR}/1.names.2 |\
    gzip >${_polap_var_oga_seeds}/3.fq.gz"

	_polap_log1 "  extracts contig sequeces from the assembly: ${_polap_var_ga_contigger_edges_fasta}"
	_polap_log2 "    input1: ${_polap_var_ga_contigger_edges_fasta}"
	_polap_log2 "    input2: ${_polap_var_mtcontigname}"
	_polap_log2 "    output: ${MTDIR}/contig.fa"
	_polap_log3_pipe "seqkit grep \
    --threads ${_arg_threads} \
    -f ${_polap_var_mtcontigname} \
		${_polap_var_ga_contigger_edges_fasta} \
		-o ${MTDIR}/contig.fa \
		2>${_polap_output_dest}"

	# 1.fq.gz: reads mapped on a single read
	# 2.fq.gz: a reduced data of 1.fq.gz
	#          or a reduced data of 1.fq.gz plus the inter-mapped reads, if any
	# single.names.2: read names from the reduced data
	#
	# v0.2.6
	# input1: single.names: single-mapped read names
	# input2: edge-number1_edge_number2.name: pair-mapped read names
	# 1.names: unique read names from single.names
	# 1.fq.gz: reads using 1.names
	# 2.fq.gz: reduced from 1.fq.gz but later replaced by some with additional pair-mapped reads
	# single.names.2: names from 2.fq.gz with only single-mapped reads or reduced names
	# 1.names.2: single.names.2 + all pair-mapped read names
	#
	# single.names.2: single-mapped reads -> read data
	# 1.names.2: reads used for organelle-genome
	# use:
	# seqtk subseq ${_polap_var_outdir_nk_fq_gz} ${MTSEEDSDIR}/1.names.2 | gzip >${MTSEEDSDIR}/3.fq.gz

	_polap_log0 "assembling the organelle-genome using single-mapped or pair-mapped reads ..."

	_polap_utility_get_contig_length \
		"${MTDIR}/contig.fa" \
		"${MTDIR}/contig_total_length.txt"
	local CONTIG_LENGTH=$(<"${MTDIR}/contig_total_length.txt")
	_polap_log2_cat "${_polap_var_oga_contig}/contig_total_length.txt"
	# CONTIG_LENGTH=$(seqkit stats -Ta "${MTDIR}"/contig.fa | csvtk cut -t -f "sum_len" | csvtk del-header)
	# echo "$CONTIG_LENGTH" >"${MTDIR}"/contig_total_length.txt

	_polap_log1 "  organelle genome size based on the seed contig selection: $CONTIG_LENGTH"

	for _i_seed in 1 2 3; do

		if [[ "${_i_seed}" = "2" ]]; then
			_polap_log0 "assembling the organelle-genome using single-mapped reads ..."
		elif [[ "${_i_seed}" = "3" ]]; then
			_polap_log0 "assembling the organelle-genome using pair-mapped reads ..."
		elif [[ "${_i_seed}" = "1" ]]; then
			_polap_log0 "assembling the organelle-genome using all single-mapped reads ..."
		fi

		_polap_log3_cmd mkdir -p "${_polap_var_oga_seeds}/${_i_seed}"

		if [[ "${_arg_menu[1]}" = "polishing" ]]; then
			_polap_log1 "  executing the organelle-genome long-read flye polishing assembly on ${_arg_jnum} ..."
			_polap_log1 "    input1: ${_polap_var_oga_seeds}/${_i_seed}.fq.gz"
			_polap_log1 "    output: ${_polap_var_oga_seeds}/${_i_seed}/assembly_graph.gfa"
		else
			_polap_log1 "  executing the organelle-genome assembly using flye on ${_arg_jnum} ..."
			_polap_log1 "    input1: ${_polap_var_oga_seeds}/${_i_seed}.fq.gz"
			_polap_log1 "    output: ${_polap_var_oga_seeds}/${_i_seed}/30-contigger/graph_final.gfa"
		fi

		local _command1="flye \
      ${_arg_flye_data_type} \
      ${_polap_var_oga_seeds}/${_i_seed}.fq.gz \
		  --out-dir ${_polap_var_oga_seeds}/${_i_seed} \
		  --threads ${_arg_threads} \
		  --asm-coverage ${_arg_flye_asm_coverage} \
		  --genome-size $CONTIG_LENGTH"
		if [[ "${_arg_menu[1]}" = "polishing" ]]; then
			_command1+=" \
		  --resume"
		else
			_command1+=" \
		  --stop-after contigger"
		fi
		_command1+=" \
		  2>${_polap_output_dest}"
		_polap_log3_pipe "${_command1}"
	done

	_polap_log1 NEXT: $(basename "$0") flye2 -o "${_arg_outdir}" -j "${_arg_jnum}" -t "${_arg_threads}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Executes Flye for an organelle-genome assembly
#
# Arguments:
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   -t ${_arg_threads}: the number of CPU cores
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${MTDIR}/contig.fa
#   ${MTSEEDSDIR}/2.fq.gz
# Outputs:
#   ${MTDIR}/contig_total_length.txt
#   ${MTDIR}/30-contigger/contigs.fasta
#   ${MTDIR}R}/30-contigger/contigs_stats.txt
#   ${MTDIR}R}/30-contigger/graph_final.fasta
#   ${MTDIR}/30-contigger/graph_final.gfa
################################################################################
function _run_polap_x-v0.2.6-flye2() { # executes Flye for an organelle-genome assembly
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	local MTDIR="${_arg_outdir}"/${_arg_jnum}
	MTSEEDSDIR="${MTDIR}"/seeds

	help_message=$(
		cat <<HEREDOC
# Executes Flye for an organelle-genome assembly
# Arguments:
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   -t ${_arg_threads}: the number of CPU cores
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${MTDIR}/contig.fa
#   ${MTSEEDSDIR}/2.fq.gz
# Outputs:
#   ${MTDIR}/contig_total_length.txt
#   ${MTDIR}/30-contigger/contigs.fasta
#   ${MTDIR}/30-contigger/contigs_stats.txt
#   ${MTDIR}/30-contigger/graph_final.fasta
#   ${MTDIR}/30-contigger/graph_final.gfa
Example: $(basename "$0") ${_arg_menu[0]} -j <arg>
Example: $(basename "$0") ${_arg_menu[0]} -j <arg> polishing
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	if [[ "${_arg_menu[1]}" = "polishing" ]]; then
		_polap_log0 "flye2 assembly with polishing on the reads mapped on seed contigs ..."
	else
		_polap_log0 "flye2 assembly on the reads mapped on seed contigs ..."
	fi
	_polap_log1 "  input1: ${MTDIR}/contig.fa"
	_polap_log1 "  input2: ${MTSEEDSDIR}/3.fq.gz"

	if [ ! -s "${MTDIR}/contig.fa" ]; then
		echoall "ERROR: no selected-contig file [${MTDIR}/contig.fa]"
		echoerr "SUGGESTION: select-reads"
		exit $EXIT_SUCCESS
	fi

	if [ ! -s "${MTSEEDSDIR}/3.fq.gz" ]; then
		echoall "ERROR: no long-read file: ${MTSEEDSDIR}/3.fq.gz"
		echoerr "SUGGESTION: select-reads"
		exit $EXIT_SUCCESS
	fi

	_polap_utility_get_contig_length \
		"${_polap_var_oga}/contig.fa" \
		"${_polap_var_oga}/contig_total_length.txt"
	local CONTIG_LENGTH=$(<"${_polap_var_oga}/contig_total_length.txt")
	_polap_log2_cat "${_polap_var_oga}/contig_total_length.txt"

	# CONTIG_LENGTH=$(seqkit stats -Ta "${MTDIR}"/contig.fa | csvtk cut -t -f "sum_len" | csvtk del-header)
	# echo "$CONTIG_LENGTH" >"${MTDIR}"/contig_total_length.txt
	_polap_log1 "  organelle genome size based on the seed contig selection: $CONTIG_LENGTH"

	if [[ "${_arg_menu[1]}" = "polishing" ]]; then
		_polap_log1 "  executing the organelle-genome long-read flye polishing assembly on ${_arg_jnum} ..."
		_polap_log1 "    input1: ${MTSEEDSDIR}/3.fq.gz"
		_polap_log1 "    output: ${MTDIR}/assembly_graph.gfa"
	else
		_polap_log1 "  executing the organelle-genome assembly using flye on ${_arg_jnum} ..."
		_polap_log1 "    input1: ${MTSEEDSDIR}/3.fq.gz"
		_polap_log1 "    output: ${MTDIR}/30-contigger/graph_final.gfa"
	fi
	local _command1="flye \
    --nano-raw ${MTSEEDSDIR}/3.fq.gz \
		--out-dir ${MTDIR} \
		--threads ${_arg_threads} \
		--asm-coverage ${_arg_flye_asm_coverage} \
		--genome-size $CONTIG_LENGTH"
	if [[ "${_arg_menu[1]}" = "polishing" ]]; then
		_command1+=" \
		  --resume"
	else
		_command1+=" \
		  --stop-after contigger"
	fi
	_command1+=" \
		2>${_polap_output_dest}"
	_polap_log3_pipe "${_command1}"

	if [[ "${_arg_menu[1]}" = "polishing" ]]; then
		_polap_log0 "  output: the long-read polished assembly graph: $PWD/${_polap_var_oga_assembly_graph_gfa}"
		_polap_log0 "  DO: extract a draft organelle genome sequence (mt.0.fasta)"
		_polap_log0 "      from the polished assembly graph before short-read polishing"
		_polap_log1 "NEXT: $(basename "$0") prepare-polishing -a ${_arg_short_read1} -b ${_arg_short_read2}"
	else
		_polap_log0 "  output: the assembly fasta: ${_polap_var_oga_contigger_edges_fasta}"
		_polap_log0 "  output: the assembly graph: ${_polap_var_oga_contigger_edges_gfa}"
		jnum_next=$((_arg_jnum + 1))
		_polap_log1 "  create and edit ${_arg_outdir}/${_arg_jnum}/mt.contig.name-${jnum_next}"
		_polap_log1 "NEXT: $(basename "$0") assemble2 -o ${_arg_outdir} -i ${_arg_jnum} -j ${jnum_next}"
		_polap_log1 or you could finish with Flye organelle-genome assembly with its polishing stage.
		_polap_log1 "NEXT: $(basename "$0") flye-polishing -o ${_arg_outdir} -j $${_arg_jnum}"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
