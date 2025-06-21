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
#
# Note: deleted these or move to run-polap-function-oga-v0.2.6.sh
# function _run_polap_x-v0.3.7-collect-reads() { # replaced by select-reads
# function _run_polap_x-v0.2.6-select-reads() { # selects reads mapped on a genome assembly
# function _run_polap_x-select-reads { # selects reads mapped on a genome assembly
# function _run_polap_x-v0.2.6-flye2() { # executes Flye for an organelle-genome assembly
#
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
# create a dummy assembly for preparing seed contigs
################################################################################
function _run_polap_prepare-seeds { # prepare seed contigs in a not usual way
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
# Prepare contig seeds from multiple assemblies or sequences
#
# Arguments:
#   -i number1,number2,...
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#
#   -f contig sequence file
# Inputs:
#   mulitple assemblies
#   sequence file
# Outputs:
#   a dummy assembly
Example: $(basename "$0") ${_arg_menu[0]} -i 4,6 -j 7
Example: $(basename "$0") ${_arg_menu[0]} -f seed_contigs.fa
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	_polap_log0 "preparing seed contigs ..."
	_polap_log1 "  assembly: ${_arg_inum} (source) -> ${_arg_jnum} (target) ..."
	# if [[ -d "${_arg_outdir}/${_arg_jnum}" ]]; then
	# 	_polap_log0 "ERROR: target assembly number ${_arg_jnum} already exists."
	# 	return $RETURN_SUCCESS
	# fi
	mkdir -p "${_arg_outdir}/${_arg_jnum}"
	if [[ "${_arg_final_assembly}" == "mt.1.fa" ]]; then
		_polap_log0 "  -i must have two or more assembly numbers"

		# Initialize a counter
		local counter=1

		# Output file for concatenated FASTA
		local output_file="${_arg_outdir}/${_arg_jnum}/seed_contigs.fa"

		# Clear or create the output file
		>"$output_file"
		# Set IFS to a comma to split the string by commas
		IFS=','

		# Use a for loop to iterate over each number
		for number in ${_arg_inum}; do
			local fasta_file="${_arg_outdir}/${number}/assembly.fasta"

			if [[ ! -s "${fasta_file}" ]]; then
				_polap_log0 "ERROR: no such file: ${fasta_file}"
				return $RETURN_SUCCESS
			fi

			# Process each file, renaming sequence IDs
			awk -v counter="$counter" '
    /^>/ {
        # Print a new sequence ID with "edge_" followed by the current counter
        print ">edge_" counter
        counter++
        next
    }
    # Print the sequence lines as they are
    {
        print
    }' "$fasta_file" >>"$output_file"

			# Update the counter based on the number of sequences processed in the current file
			sequences_in_file=$(grep -c "^>" "$fasta_file")
			counter=$((counter + sequences_in_file))
		done

		# Reset IFS to default (optional, but good practice)
		unset IFS

		_polap_log0 "contig sequence file: $output_file"

	else
		# cp "${_arg_final_assembly}" "${_arg_outdir}/${_arg_jnum}/seed_contigs.fa"

		# Input FASTA file
		input_file="${_arg_final_assembly}"
		# Output FASTA file
		output_file="${_arg_outdir}/${_arg_jnum}/seed_contigs.fa"

		# Initialize a counter
		local counter=1

		# Process the input file and rename sequence IDs
		awk -v counter="$counter" '
/^>/ {
    # Print a new sequence ID with "edge_" followed by the current counter
    print ">edge_" counter
    counter++
    next
}
# Print the sequence lines as they are
{
    print
}' "$input_file" >"$output_file"

		_polap_log0 "  ${output_file} is the contig sequence file."
	fi

	_polap_log3_cmd mkdir "${_polap_var_oga_contigger}"
	_polap_log3_cmd cp "${output_file}" "${_polap_var_oga_contigger_edges_fasta}"
	local _arg_knum=$((_arg_jnum + 1))
	seqkit seq -n -i "${output_file}" >"${_polap_var_oga}/mt.contig.name-${_arg_knum}"
	_polap_log0_cat "${_polap_var_oga}/mt.contig.name-${_arg_knum}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Map reads on a genome assembly using minimap2.
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
#   _polap_var_mtcontigname="${_polap_var_ga}/mt.contig.name-${_arg_jnum}"
#   ${_polap_var_mtcontigname} <- ${_arg_inum} and ${_arg_jnum}
#   ${_polap_var_mtcontigname} <- ${_arg_inum} and ${_arg_jnum}
#   ${_polap_var_ga_contigger_edges_fasta} <- ${_arg_inum}
#   ${_polap_var_outdir_nk_fq_gz}
#   ${_polap_var_outdir_lk_fq_gz}
#   ${_polap_var_outdir_l_fq_gz}
#
# Outputs:
#   _polap_var_oga_reads="${_polap_var_oga}/02-reads"
#   ${_polap_var_oga_reads}/contig.fa
#   ${_polap_var_oga_reads}/contig.paf
#   ${_polap_var_oga_reads}/contig.tab
################################################################################
function _run_polap_map-reads { # selects reads mapped on a genome assembly
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
# Map long reads on a Flye genome assembly.
#
# Arguments:
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number (or 0)
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   -l ${_arg_long_reads}: long-read data default:${_polap_var_outdir_lk_fq_gz}
# Inputs:
#   ${_polap_var_mtcontigname}
#   ${_polap_var_ga_contigger_edges_fasta}
#   ${_polap_var_outdir_lk_fq_gz} or ${_polap_var_outdir_nk_fq_gz}
# Outputs:
#   ${_polap_var_oga_contig}/contig.fa
#   ${_polap_var_oga_contig}/contig.paf
#   ${_polap_var_oga_contig}/contig.tab
Example: $(basename "$0") ${_arg_menu[0]} -i ${_arg_inum} -j ${_arg_jnum}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	_polap_log0 "mapping long-read data on the seed contigs ..."
	_polap_log1 "  assembly: ${_arg_inum} (source) -> ${_arg_jnum} (target) ..."
	_polap_log1 "  input1: ${_polap_var_mtcontigname}"
	_polap_log1 "  input2: ${_polap_var_ga_contigger_edges_fasta}"

	if [[ -s "${_polap_var_oga_contig}/contig.tab" ]] && [[ "${_arg_redo}" == "on" ]]; then
		_polap_log0 "  found: ${_polap_var_oga_reads}/contig.tab, so skipping mapping long-read data ..."
		return
	fi

	if [ ! -s "${_polap_var_mtcontigname}" ]; then
		_polap_log0 "ERROR: no such mt.contig.name file: ${_polap_var_mtcontigname}"
		exit $EXIT_ERROR
	fi

	if [ ! -s "${_polap_var_ga_contigger_edges_fasta}" ]; then
		_polap_log0 "ERROR: no assembly fasta file: ${_polap_var_ga_contigger_edges_fasta}"
		exit $EXIT_ERROR
	fi

	if [ -d "${_polap_var_oga}" ] && [ "${_arg_redo}" = "off" ]; then
		while true; do
			read -r -p "Folder [${_polap_var_oga}] already exists. Do you want to replace it? [y/n] " yn
			case $yn in
			[Yy]*)
				_polap_log3_cmd rm -rf "${_polap_var_oga}"
				_polap_log0 "  folder ${_polap_var_oga} is deleted."
				break
				;;
			[Nn]*)
				_polap_log0 "  folder ${_polap_var_oga} is not replaced."
				_polap_log0 "  you might want a new mt.contig.name file for Flye2 step."
				exit $EXIT_FAIL
				;;
			*) _polap_log0 "Please answer yes or no." ;;
			esac
		done
	else
		_polap_log3_cmd rm -rf "${_polap_var_oga}"
		_polap_log1 "  ${_polap_var_oga} is deleted if there is one."
	fi

	_polap_log1 "  creates ${_polap_var_oga_contig}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_contig}"

	_polap_log1 "  determines which long-read data to use ..."

	local _source_long_reads_fq=""
	_polap_oga_determine-long-read-file _source_long_reads_fq

	if ! _polap_gfatools-gfa2fasta; then
		_polap_error_message $?
		return ${_POLAP_ERR_MENU_MAP_READS}
	fi

	# _polap_log3_cmd ln -s "$PWD/${_polap_var_outdir_nk_fq_gz}" -t "${_polap_var_oga}"
	# _polap_log3_cmd ln -s "$PWD/${_polap_var_outdir_lk_fq_gz}" -t "${_polap_var_oga}"
	_polap_log1 "  extracts contig sequeces from the assembly: ${_polap_var_ga_contigger_edges_fasta}"
	_polap_log2 "    input1: ${_polap_var_ga_contigger_edges_fasta}"
	_polap_log2 "    input2: ${_polap_var_mtcontigname}"
	_polap_log2 "    output: ${_polap_var_oga_contig}/contig.fa"
	_polap_log3_pipe "seqkit grep \
    -f ${_polap_var_mtcontigname} \
		${_polap_var_ga_contigger_edges_fasta} \
		-o ${_polap_var_oga_contig}/contig.fa \
		2>${_polap_output_dest}"

	# we could circularize a single contig.
	local contig_count=$(wc -l <"${_polap_var_mtcontigname}")
	_polap_log1 "  number of seed contigs in ${_polap_var_mtcontigname}: ${contig_count}"
	if [[ ${_arg_circularize} == "on" ]]; then
		if [ "$contig_count" -eq 1 ]; then

			_polap_log1 "  circularizes the single seed contig ..."
			_polap_log2 "    input1: ${_polap_var_oga_contig}"
			_polap_log2 "    input2: ${_polap_var_oga_contig}/contig.fa"
			_polap_log2 "    input3: ${_polap_var_mtcontigname}"
			_polap_log2 "    output1: ${_polap_var_mtcontigname}-backup"
			_polap_log2 "    output2: (new) ${_polap_var_mtcontigname}"
			_polap_log2 "    output3: ${_polap_var_oga_contig}/contig.fa-backup"
			_polap_log2 "    output4: (new) ${_polap_var_oga_contig}/contig.fa"
			_polap_log3_cmd bash "${_POLAPLIB_DIR}/run-polap-sh-half-cut.sh" "${_polap_var_oga_contig}" "${_polap_var_oga_contig}/contig.fa" "${_polap_var_mtcontigname}"

			# cp "${_polap_var_oga}/contig.fa" "${_polap_var_oga}/contig.fa-backup"
			# seqkit fx2tab --length --name "${_polap_var_oga}"/contig.fa -o "${_polap_var_oga}"/contig.fa.len >/dev/null 2>&1
			# A=$(cut -f2 "${_polap_var_oga}"/contig.fa.len)
			# B=$(echo "scale=0; $A/2" | bc)
			# C=$((B + 1))
			# seqkit subseq -r 1:"$B" "${_polap_var_oga}"/contig.fa -o "${_polap_var_oga}"/c1.fa >/dev/null 2>&1
			# seqkit subseq -r "$C":"$A" "${_polap_var_oga}"/contig.fa -o "${_polap_var_oga}"/c2.fa >/dev/null 2>&1
			# cat "${_polap_var_oga}"/c?.fa | seqkit replace -p '.+' -r 'edge_{nr}' -o "${_polap_var_oga}"/contig.fa >/dev/null 2>&1
			# cp "${_polap_var_mtcontigname}" "${_polap_var_mtcontigname}"-backup
			# echo -e "edge_1\nedge_2" >"${_polap_var_mtcontigname}"

			_polap_log1 "    creating new ${_polap_var_oga_reads}/contig.fa and ${_polap_var_mtcontigname}"

		else
			_polap_log0 "Not implemented yet!"
			exit $EXIT_ERROR
		fi
	fi

	_polap_log1 "  finding the length of all seed contigs"
	_polap_utility_get_contig_length \
		"${_polap_var_oga_contig}/contig.fa" \
		"${_polap_var_oga_contig}/contig_total_length.txt"
	local CONTIG_LENGTH=$(<"${_polap_var_oga_contig}/contig_total_length.txt")
	_polap_log2_cat "${_polap_var_oga_contig}/contig_total_length.txt"
	# local CONTIG_LENGTH=$(seqkit stats -Ta "${_polap_var_oga_contig}"/contig.fa | csvtk cut -t -f "sum_len" | csvtk del-header)
	# echo "$CONTIG_LENGTH" >"${_polap_var_oga_contig}"/contig_total_length.txt
	local _contig_length_bp=$(_polap_utility_convert_bp ${CONTIG_LENGTH})
	_polap_log1 "    organelle genome size based on the seed contig selection: ${_contig_length_bp}"

	_polap_log1 "  mapping long-read data on the seed contigs using minimap2 ..."
	_polap_log2 "    input1: ${_polap_var_oga_contig}/contig.fa"
	_polap_log2 "    input2: ${_source_long_reads_fq}"
	_polap_log2 "    output: ${_polap_var_oga_contig}/contig.paf"
	if [[ -s "${_polap_var_oga_contig}"/contig.paf ]] && [[ "${_arg_redo}" = "off" ]]; then
		_polap_log1 "  found: ${_polap_var_oga_reads}/contig.paf, skipping the minimap2 mapping step ..."
	else
		_polap_log3_pipe "minimap2 -cx \
      ${_arg_minimap2_data_type} \
      ${_polap_var_oga_contig}/contig.fa \
      ${_source_long_reads_fq} \
      -t ${_arg_threads} \
      -o ${_polap_var_oga_contig}/contig.paf \
      >${_polap_output_dest} 2>&1"
	fi

	_polap_log1 "  converting PAF to TAB ..."
	_polap_log2 "    input1: ${_polap_var_oga_contig}/contig.paf"
	_polap_log2 "    output: ${_polap_var_oga_contig}/contig.tab"
	# cut -f1-11 "${_polap_var_oga}"/contig.paf | awk -v minlength="${_arg_min_read_length}" '{if ($2>=minlength) {print}}' >"${_polap_var_oga}"/contig.tab
	_polap_log3_cmd bash "${_POLAPLIB_DIR}/run-polap-sh-minimap2-paf2tab.sh" "${_arg_min_read_length}" "${_polap_var_oga_contig}/contig.paf" "${_polap_var_oga_contig}/contig.tab"

	_polap_log1 "NEXT: $(basename "$0") reads -o ${_arg_outdir} -i ${_arg_inum} -j ${_arg_jnum}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
	return

	# put the backup to the original
	if [[ ${_arg_circularize} == "on" ]]; then
		_polap_log1 "  circularize option: putting seed contigs back to the original mt.contig.name and ${_polap_var_oga}/contig.fa"
		if [[ -s "${_polap_var_mtcontigname}"-backup ]]; then
			mv "${_polap_var_mtcontigname}"-backup "${_polap_var_mtcontigname}"
			mv "${_polap_var_oga}/contig.fa-backup" "${_polap_var_oga}/contig.fa"
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
# Test organelle assemblies on a range of w values.
#
# Arguments:
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   -s, --select-read-range <start,end,count>
#   --start-index <index>: to start at somewhere not start
# Inputs:
#   ${_polap_var_mtcontigname}
#   ${_polap_var_ga_contigger_edges_fasta}
#   input long read data to use (priority order):
#     1. ${_polap_var_outdir_lk_fq_gz}
#     2. ${_polap_var_outdir_nk_fq_gz}
#     3. ${_arg_long_reads}
# Outputs:
#   ${_polap_var_oga_contig}: map-reads output
#   ${_polap_var_oga_reads}: mapped read names
#   ${_polap_var_oga_seeds}: mapped reads
#   ${_polap_var_oga_subsample}: subsample of the mapped reads
#   ${_polap_var_oga_flye}: assemblies
#   ${_polap_var_oga_summary}: summary
#   ${_polap_var_oga_plot}: plot or table using range in contig folder and summary
#
# Menus:
#   ptgaul-flye-asm-coverage -> ptgaul --select-read-range 30,90,3
#   ptgaul-intra-base-ratio  -> ptgaul --select-read-range 0.1,0.9,5
#   ptgaul-intra-base-length -> ptgaul --select-read-range 3000,27000,5
#   single-intra-read-ratio   -> single
#   single-intra-base-ratio   -> single
#   single-intra-base-length  -> single
#   combined-intra-read-ratio   -> combined
#   combined-intra-base-ratio   -> combined
#   combined-intra-base-length  -> combined
#   combined-inter-base-ratio   -> combined
#   combined-inter-base-length  -> combined
#   bridge-inter-base-length  -> bridge
#   combined-bridge-read-length -> combined
# View:
#   ptgaul-intra-base-length
#   single-intra-base-length
#   polap-rw-base-length
################################################################################
function _run_polap_test-reads { # selects reads mapped on a genome assembly
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
# Test organelle assemblies on a range of w values.
#
# Arguments:
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   -s, --select-read-range <start,end,count>
#   --start-index <index>: to start at somewhere not start
# Inputs:
#   ${_polap_var_mtcontigname}
#   ${_polap_var_ga_contigger_edges_fasta}
#   input long read data to use (priority order): 
#     1. ${_polap_var_outdir_lk_fq_gz}
#     2. ${_polap_var_outdir_nk_fq_gz}
#     3. ${_arg_long_reads}
# Outputs:
#   ${_polap_var_oga_contig}: map-reads output
#   ${_polap_var_oga_reads}: mapped read names
#   ${_polap_var_oga_seeds}: mapped reads
#   ${_polap_var_oga_subsample}: subsample of the mapped reads
#   ${_polap_var_oga_flye}: assemblies
#   ${_polap_var_oga_summary}: summary
#   ${_polap_var_oga_plot}: plot or table using range in contig folder and summary
#
# Menus:
#   ptgaul-intra-base-length -> ptgaul --select-read-range 3000,27000,5
#   single-intra-base-length  -> single
#   polap-rw-base-length -> combined
#   bridge-inter-base-length  -> bridge
# View:
#   ptgaul-intra-base-length
#   single-intra-base-length
#   polap-rw-base-length
Example: $(basename "$0") ${_arg_menu[0]} [ptgaul-intra-base-length] --select-read-range 3000,27000,5
Example: $(basename "$0") ${_arg_menu[0]} single-intra-base-length -i 1 -j 2
Example: $(basename "$0") ${_arg_menu[0]} polap-rw-base-length --select-read-range 3000,27000,5
Example: $(basename "$0") ${_arg_menu[0]} bridge-inter-base-length
Example: $(basename "$0") ${_arg_menu[0]} polap-rw-base-length --select-read-range 3000,27000,5 --start-index 3
Example: $(basename "$0") ${_arg_menu[0]} view polap-rw-base-length -i 2 -j 3
Example: $(basename "$0") ${_arg_menu[0]} report ptgaul --report-x 3000,5000,7000,9000,11000
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	if [[ "${_arg_menu[1]}" == "preview" ]]; then
		seqkit stats -Ta "${_polap_var_oga_contig}/contig.fa" >&3

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	_polap_log0 "test read-selection ..."

	declare -A _menu_type_dict
	_menu_type_dict["ptgaul"]="ptgaul-intra-base-length"
	_menu_type_dict["intra"]="single-intra-base-length"
	_menu_type_dict["polap"]="polap-rw-base-length"
	_menu_type_dict["bridge"]="bridge-inter-base-length"

	declare -A _polap_opt_dict
	_polap_opt_dict["ptgaul-flye-asm-coverage"]="ptgaul"
	_polap_opt_dict["ptgaul-intra-base-ratio"]="ptgaul"
	_polap_opt_dict["ptgaul-intra-base-length"]="ptgaul"
	_polap_opt_dict["single-intra-base-ratio"]="single"
	_polap_opt_dict["single-intra-base-ratio"]="single"
	_polap_opt_dict["single-intra-base-length"]="single"
	_polap_opt_dict["combined-intra-read-ratio"]="combined"
	_polap_opt_dict["combined-intra-base-ratio"]="combined"
	_polap_opt_dict["combined-intra-base-length"]="combined"
	_polap_opt_dict["combined-inter-base-ratio"]="combined"
	_polap_opt_dict["combined-inter-base-length"]="combined"
	_polap_opt_dict["bridge-inter-base-length"]="pair"
	_polap_opt_dict["polap-rw-base-length"]="combined"
	_polap_opt_dict["combined-bridge-length"]="combined"
	_polap_opt_dict["ptgaul-reads"]="ptgaul"
	_polap_opt_dict["intra-reads"]="single"
	_polap_opt_dict["bridge-reads"]="pair"
	_polap_opt_dict["polap-reads"]="combined"

	if [[ "${_arg_menu[1]}" == "view" ]]; then
		# use the summary to plot the number of bases
		# use the summary to plot the number of fragments
		if [[ "${_arg_menu[2]}" == "outfile" ]]; then
			_arg_menu[2]="ptgaul-intra-base-length"
			_polap_log0 "  default set to read-selection: ${_arg_menu[2]}"
		fi
		local _pread_sel="${_arg_menu[2]}"
		if [[ -v _menu_type_dict["${_pread_sel}"] ]]; then
			_pread_sel="${_menu_type_dict["${_pread_sel}"]}"
		fi

		if [[ -v _polap_opt_dict["${_pread_sel}"] ]]; then
			local _read_names=${_polap_opt_dict["${_pread_sel}"]}
			_polap_log0 "  read names file: ${_read_names}"
		else
			_polap_log0 "ERROR: no such name file for ${_pread_sel}"
			return
		fi

		last_number=$(ls ${_polap_var_oga_summary}/${_pread_sel}/*.bases | sed 's#.*/##' | sed 's/\.bases//' | sort -n | tail -1)
		for i in $(seq 0 "$last_number"); do
			cat "${_polap_var_oga_summary}/${_pread_sel}/$i.bases"
		done >"${_polap_var_oga_summary}/${_pread_sel}/bases.txt"

		for i in $(seq 0 "$last_number"); do
			cat "${_polap_var_oga_summary}/${_pread_sel}/$i.fragments"
		done >"${_polap_var_oga_summary}/${_pread_sel}/fragments.txt"

		for i in $(seq 0 "$last_number"); do
			cat "${_polap_var_oga_summary}/${_pread_sel}/$i.depth"
		done >"${_polap_var_oga_summary}/${_pread_sel}/depth.txt"

		# Define filenames
		local file1="${_polap_var_oga_contig}/${_pread_sel}.txt"
		local file2="${_polap_var_oga_summary}/${_pread_sel}/bases.txt"
		local file3="${_polap_var_oga_summary}/${_pread_sel}/fragments.txt"
		local file4="${_polap_var_oga_summary}/${_pread_sel}/depth.txt"

		if [[ ! -s "$file2" ]]; then
			_polap_log0 "ERROR: no such file: $file2"
			return
		fi
		local output_file="${_polap_var_oga_summary}/${_pread_sel}/summary.tsv"

		# Create header for output file with an index column
		echo -e "index\tsize\tbases\tfragments\tdepth" >"$output_file"

		# Split values from file1 into lines, combine with file2 and file3, and output to file
		awk '{ for (i=1; i<=NF; i++) print $i }' "$file1" | paste - "$file2" "$file3" "$file4" | awk '{ print NR-1 "\t" $0 }' >>"$output_file"

		_polap_log0_column "$output_file"

		_polap_log1 "  selecting read names for ${_pread_sel}"
		_polap_log2 "    input1: ${output_file}"
		_polap_log2 "    output: ${_polap_var_oga_plot}/${_pread_sel}/summary.pdf"
		_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/run-polap-r-test-reads-bar-graph.R \
			-i ${output_file} \
			--out ${_polap_var_oga_plot}/${_pread_sel}-summary.pdf \
			>${_polap_output_dest} 2>&1"

		_polap_log0 "See: ${_polap_var_oga_plot}/${_pread_sel}-summary.pdf"

		local _project_dir=$(grep 'Project:' ${LOG_FILE} | tail -1 | cut -d: -f4 | head -n 1 | xargs)
		local _species=$(grep 'species:' ${LOG_FILE} | tail -1 | cut -d: -f4 | head -n 1 | xargs)

		local _md_file="${_polap_var_oga_plot}/figure-${_pread_sel}.md"
		printf "%s\n\n" "# mtDNA of ${_species}: ${_pread_sel}" >"${_md_file}"
		read -a restored_array <"${_polap_var_oga_contig}/${_pread_sel}.txt"
		local array_length=${#restored_array[@]}
		# Iterate over the array using an index
		for ((i = ${_arg_start_index}; i < array_length; i++)); do
			local _test_value="${restored_array[i]}"
			_polap_log3_pipe "cp -f ${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph.png \
			  ${_polap_var_oga_flye}/${_pread_sel}/${_test_value}-graph_final.png 2>${_polap_output_dest}"
			if [[ "${_pread_sel}" =~ -reads$ ]]; then
				_polap_log3_pipe "cp ${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa \
			    ${_polap_var_oga_plot}/${_pread_sel}/${i}-${_test_value}-graph_final.gfa"
			else
				_polap_log3_pipe "cp ${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa \
			    ${_polap_var_oga_plot}/${_pread_sel}/${_test_value}-graph_final.gfa"
				printf "![${_test_value}](%s){ width=13%% }\n" ${_project_dir}/${_polap_var_oga_plot}/${_pread_sel}/${_test_value}-graph_final.png >>"${_md_file}"
			fi
		done
		printf "\n![Numbers of bases and fragments: ${_pread_sel} of ${_species}](%s)\n" ${_project_dir}/${_polap_var_oga_plot}/${_pread_sel}-summary.pdf >>"${_md_file}"
		printf "\n\n" >>"${_md_file}"
		echo "\newpage" >>"${_md_file}"
		printf "\n\n" >>"${_md_file}"

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	if [[ "${_arg_menu[1]}" == "report" ]]; then
		local _pread_sel="${_arg_menu[2]}"
		if [[ -v _menu_type_dict["${_pread_sel}"] ]]; then
			_pread_sel="${_menu_type_dict["${_pread_sel}"]}"
		fi
		local output_file="${_polap_var_oga_summary}/${_pread_sel}/summary.tsv"
		_polap_log1 "  selecting read names for ${_pread_sel}"
		_polap_log2 "    input1: ${output_file}"
		_polap_log2 "    output: ${_polap_var_oga_plot}/${_pread_sel}/summary.pdf"
		if [[ "${_arg_report_x_is}" == "off" ]]; then
			_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/run-polap-r-test-reads-bar-graph.R \
			-i ${output_file} \
			--out ${_polap_var_oga_plot}/${_pread_sel}-summary.pdf \
			>${_polap_output_dest} 2>&1"
		else
			_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/run-polap-r-test-reads-bar-graph.R \
			-i ${output_file} \
			-s ${_arg_report_x} \
			--out ${_polap_var_oga_plot}/${_pread_sel}-summary.pdf \
			>${_polap_output_dest} 2>&1"
		fi

		_polap_log0 "See: ${_polap_var_oga_plot}/${_pread_sel}-summary.pdf"

		# Convert the string to an array by changing the IFS to a comma
		IFS=',' read -r -a restored_array <<<"${_arg_report_x}"

		local _project_dir="project-xyz"
		local _md_file="${_polap_var_oga_plot}/figure-${_pread_sel}.md"
		printf "%s\n\n" "# mtDNA ${_pread_sel}" >"${_md_file}"
		# read -a restored_array <"${_polap_var_oga_contig}/${_pread_sel}.txt"
		local array_length=${#restored_array[@]}
		# Iterate over the array using an index
		for ((i = ${_arg_start_index}; i < array_length; i++)); do
			local _test_value="${restored_array[i]}"
			printf "![${_test_value}](%s){ width=13%% }\n" ${_project_dir}/${_polap_var_oga_plot}/${_pread_sel}/${_test_value}-graph_final.png >>"${_md_file}"
		done
		printf "\n![Numbers of bases and fragments: ${_pread_sel}](%s)\n" ${_project_dir}/${_polap_var_oga_plot}/${_pread_sel}-summary.pdf >>"${_md_file}"
		printf "\n\n" >>"${_md_file}"
		echo "\newpage" >>"${_md_file}"
		printf "\n\n" >>"${_md_file}"
		_polap_log0 "See: ${_md_file}"

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	##############################################################################
	# Set the default submenu option for the main menu
	#   _pread_sel: submenu
	#   _read_names: read name file to use (ptgaul, single, combined)
	if [[ "${_arg_menu[1]}" == "infile" ]]; then
		_arg_menu[1]="ptgaul-intra-base-length"
		_polap_log0 "  default set to read-selection: ${_arg_menu[1]}"
	fi
	local _pread_sel="${_arg_menu[1]}"
	if [[ -v _menu_type_dict["${_pread_sel}"] ]]; then
		_pread_sel="${_menu_type_dict["${_pread_sel}"]}"
	fi

	if [[ -v _polap_opt_dict["${_pread_sel}"] ]]; then
		local _read_names=${_polap_opt_dict["${_pread_sel}"]}
		_polap_log0 "  read names file: ${_read_names}"
	else
		die "No such name file for ${_pread_sel}"
	fi

	_polap_log1 "  argument1: ${_pread_sel}"
	_polap_log1 "  argument2: ${_read_names}"

	##############################################################################
	# Main input files: seed contigs and long-read data to map on the contigs
	#   The Minimap2 mapping should have been completed
	_polap_log1 "  input1: ${_polap_var_oga_contig}"

	check_folder_existence "${_polap_var_oga_contig}"

	_polap_utility_get_contig_length \
		"${_polap_var_oga_contig}/contig.fa" \
		"${_polap_var_oga_contig}/contig_total_length.txt"
	local CONTIG_LENGTH=$(<"${_polap_var_oga_contig}/contig_total_length.txt")
	_polap_log2_cat "${_polap_var_oga_contig}/contig_total_length.txt"

	_polap_log1 "  determines which long-read data to use ..."
	local _source_long_reads_fq=""
	_polap_oga_determine-long-read-file _source_long_reads_fq

	_polap_log1 "  input2: ${_source_long_reads_fq}"
	if [[ -s "${_source_long_reads_fq}" ]]; then
		_polap_log2 "  input1: ${_source_long_reads_fq} for a source of long-read data"
	else
		_polap_log0 "ERROR: no such file: ${_source_long_reads_fq}"
		return
	fi

	##############################################################################
	# create range files
	_polap_log1 "  input3: range of read-selection: ${_pread_sel}"
	case "${_pread_sel}" in
	*-ratio)
		_create_range_float "${_arg_select_read_range}" \
			"${_polap_var_oga_contig}/${_pread_sel}.txt"
		_polap_log2_cat "${_polap_var_oga_contig}/${_pread_sel}.txt"
		;;
	*-length | *-flye-asm-coverage | *-reads)
		_create_range "${_arg_select_read_range}" \
			"${_polap_var_oga_contig}/${_pread_sel}.txt"
		_polap_log2_cat "${_polap_var_oga_contig}/${_pread_sel}.txt"
		;;
	*)
		die "ERROR: no such read-selection: ${_pread_sel}"
		;;
	esac

	# Check if one wants to replace the previous run
	if [[ -d "${_polap_var_oga_reads}/${_pread_sel}" ]] &&
		[[ "${_arg_start_index}" -eq 0 ]] &&
		[[ "${_arg_yes}" == "off" ]] &&
		[[ "${_arg_redo}" == "off" ]]; then
		if confirm "Do you want to redo the test-reads of ${_pread_sel} on assembly ${_arg_jnum}"?; then
			_arg_redo="on"
		else
			_polap_log0 "You cancelled the operation."
			return 0
		fi
	fi

	if [[ "${_arg_redo}" == "on" ]] && [[ "${_arg_yes}" == "off" ]]; then
		_polap_log0 "  deleting the folder like ${_polap_var_oga_reads}/${_pread_sel}"
		rm -rf "${_polap_var_oga_reads}/${_pread_sel}"
		rm -rf "${_polap_var_oga_seeds}/${_pread_sel}"
		rm -rf "${_polap_var_oga_subsample}/${_pread_sel}"
		rm -rf "${_polap_var_oga_flye}/${_pread_sel}"
		rm -rf "${_polap_var_oga_summary}/${_pread_sel}"
		rm -rf "${_polap_var_oga_plot}/${_pread_sel}"
	fi

	_polap_log2 "  creating folders for read selection type: ${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_reads}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_seeds}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_subsample}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_flye}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_summary}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_plot}/${_pread_sel}"

	# Read the file contents into an array
	read -a restored_array <"${_polap_var_oga_contig}/${_pread_sel}.txt"
	local array_length=${#restored_array[@]}
	# Iterate over the array using an index
	for ((i = ${_arg_start_index}; i < array_length; i++)); do
		local _test_value="${restored_array[i]}"

		case "${_pread_sel}" in
		*-flye-asm-coverage)
			_arg_flye_asm_coverage="${restored_array[i]}"
			;;
		*-intra-base-ratio)
			die "Not implemented yet!"
			;;
		*-intra-read-ratio)
			die "Not implemented yet!"
			;;
		*-intra-base-length | ptgaul-reads | intra-reads)
			_arg_single_min="${restored_array[i]}"
			;;
		combined-inter-base-ratio)
			die "Not implemented yet!"
			;;
		combined-inter-base-length)
			_arg_pair_min="${restored_array[i]}"
			;;
		bridge-inter-base-length)
			_arg_pair_min="${restored_array[i]}"
			;;
		combined-bridge-read-length)
			_arg_bridge_min="${restored_array[i]}"
			;;
		polap-rw-base-length | polap-reads)
			_arg_single_min="${restored_array[i]}"
			_arg_pair_min="${restored_array[i]}"
			_arg_bridge_min="0"
			;;
		*)
			die "ERROR: ${_arg_menu[1]}"
			;;
		esac

		_polap_log3_cmd mkdir -p "${_polap_var_oga_reads}/${_pread_sel}/${i}"

		_polap_log0 "  (${i}) selecting read names for ${_pread_sel}: w = ${_test_value} (bp)"
		_polap_log2 "    single-min: ${_arg_single_min}"
		_polap_log2 "    pair-min: ${_arg_pair_min}"
		_polap_log2 "    bridge-min: ${_arg_bridge_min}"
		_polap_log2 "    input1: ${_polap_var_mtcontigname}"
		_polap_log2 "    input2: ${_polap_var_oga_contig}/contig.tab"
		_polap_log2 "    output: ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.names"
		if [[ "${_arg_bridge_same_strand}" == "off" ]]; then
			_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-pairs.R \
	    -m ${_polap_var_mtcontigname} \
		  -t ${_polap_var_oga_contig}/contig.tab \
		  --out ${_polap_var_oga_reads}/${_pread_sel}/${i} \
      -w ${_arg_single_min} \
		  -r ${_arg_pair_min} \
		  -x ${_arg_bridge_min} \
      --all \
		  >${_polap_output_dest} 2>&1"
		else
			_polap_log3_pipe "Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-bridge.R \
        --use-strand \
	    -m ${_polap_var_mtcontigname} \
		  -t ${_polap_var_oga_contig}/contig.tab \
		  --out ${_polap_var_oga_reads}/${_pread_sel}/${i} \
      -w ${_arg_single_min} \
		  -r ${_arg_pair_min} \
		  -x ${_arg_bridge_min} \
      --all \
		  >${_polap_output_dest} 2>&1"
		fi

		_polap_log1 "  selecting long reads for ${_pread_sel}"
		_polap_log2 "    input1: ${_source_long_reads_fq}"
		_polap_log2 "    input2: ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.names"
		_polap_log2 "    output: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"

		# Split the names
		# split "${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.names"

		if [[ -d "${_polap_var_outdir_lk_fq_gz}.split" ]]; then
			mkdir -p "${_polap_var_oga_seeds}/${_pread_sel}/${i}"
			ls "${_polap_var_outdir_lk_fq_gz}.split"/*.fq.gz |
				parallel seqtk subseq {} "${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.names" ">" "${_polap_var_oga_seeds}/${_pread_sel}/${i}/${_read_names}".{/}.fq

			# cat "${_polap_var_outdir_lk_fq_gz}.split"/*.fq.gz.fq |
			cat "${_polap_var_oga_seeds}/${_pread_sel}/${i}/${_read_names}".*.fq.gz.fq |
				gzip >"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"
			rm -rf "${_polap_var_oga_seeds}/${_pread_sel}/${i}"
			# rm -f "${_polap_var_oga_seeds}/${_pread_sel}/${i}/${_read_names}".*.fq.gz.fq
			# rm -f "${_polap_var_outdir_lk_fq_gz}.split"/*.fq.gz.fq
		else
			_polap_log3_pipe "seqtk subseq \
		    ${_source_long_reads_fq} \
		    ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.names |\
		    gzip >${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"
		fi

		# sample
		# we need a sample of the seeds so that it meets the coverage.
		_polap_log1 "  sampling reads ..."
		_polap_log2 "    input1: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"
		local _contig_length_bp=$(_polap_utility_convert_bp ${CONTIG_LENGTH})
		_polap_log2 "    input2 (seed contig size): ${_contig_length_bp}"
		_polap_utility_get_contig_length \
			"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz" \
			"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len"
		local _seeds_length=$(<"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len")
		_polap_log2_cat "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.len"
		local _seeds_length_bp=$(_polap_utility_convert_bp ${_seeds_length})
		_polap_log2 "    result1 (total size of reads mapped contigs): ${_seeds_length_bp}"
		local _expected_organelle_coverage=$((_seeds_length / CONTIG_LENGTH))
		_polap_log2 "    result2 (expected organelle coverage): ${_expected_organelle_coverage}x"
		if [[ "$_expected_organelle_coverage" -gt "${_arg_coverage_oga}" ]]; then
			if [[ "${_arg_coverage_check}" == "on" ]]; then
				local _rate=$(echo "scale=9; ${_arg_coverage_oga}/$_expected_organelle_coverage" | bc)
				_polap_log0 "  long-read data reduction by rate of ${_rate} <= COV[${_arg_coverage_oga}] / long-read organelle coverage[$_expected_organelle_coverage]"
				_polap_log1 "    sampling long-read data by ${_rate} ... wait ..."
				_polap_lib_random-get
				local _random_seed=${_polap_var_random_number}
				# local _random_seed=11
				_polap_log1 "    random seed for reducing long reads mapped on potential seed contigs: ${_random_seed}"
				_polap_log3_pipe "seqkit sample \
          -p ${_rate} \
          -s ${_random_seed} \
			    ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz \
          -o ${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz \
          2>${_polap_output_dest}"
				_polap_log3_pipe "echo ${_random_seed} >${_polap_var_oga_subsample}/${_pread_sel}/${i}.random.seed.${_random_seed}"
			else
				_polap_log0 "    no reduction of the long-read data because of the option --no-coverage-check: expected coverage: ${_expected_organelle_coverage}"
				_polap_log3_cmd ln -s $(realpath "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz") "${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz"
			fi
		else
			_polap_log0 "    no reduction of the long-read data because $_expected_organelle_coverage < ${_arg_coverage_oga}"
			_polap_log3_cmd ln -s $(realpath "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz") "${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz"
		fi

		_polap_log1 "  flye assembly for ${_pread_sel}"
		_polap_log2 "    input1: ${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz"
		_polap_log2 "    output: ${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa"
		# if [[ "${_arg_plastid}" == "on" ]]; then
		# 	CONTIG_LENGTH=$((CONTIG_LENGTH * 3))
		# fi
		local _command1="flye \
      ${_arg_flye_data_type} \
      ${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz \
		  --out-dir ${_polap_var_oga_flye}/${_pread_sel}/${i} \
		  --threads ${_arg_threads}"
		if [[ "${_arg_flye_asm_coverage}" -gt 0 ]]; then
			_command1+=" \
		  --asm-coverage ${_arg_flye_asm_coverage} \
		  --genome-size $CONTIG_LENGTH"
		fi
		if [[ "${_arg_menu[2]}" == "polishing" ]]; then
			_command1+=" \
		  --resume"
		else
			_command1+=" \
		  --stop-after contigger"
		fi
		_command1+=" \
		  2>${_polap_output_dest}"

		if [[ "${_arg_flye}" == "on" ]]; then
			_polap_log3_pipe "${_command1}"
		else
			_polap_log0 "No flye run in test-reads"
		fi

		# summary of the assembly
		# 1. count fragments
		# 2. count bases
		# 3. average depth
		if [[ -s "${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa" ]]; then
			_polap_log1 "  summary for ${_pread_sel}"
			_polap_log2 "    input1: ${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa"
			_polap_log2 "    output: ${_polap_var_oga_summary}/${_pread_sel}/${i}.fragments"
			_polap_log2 "    output: ${_polap_var_oga_summary}/${_pread_sel}/${i}.bases"
			_polap_log0 "  output: ${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa"

			ln -s $(realpath "${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa") ${_polap_var_oga_flye}/${_pread_sel}/${i}-graph_final.gfa

			_polap_log3_pipe "gfatools view \
        -S ${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa \
		    2>${_polap_output_dest} \
        >${_polap_var_oga_summary}/${_pread_sel}/${i}.gfa"

			# Count the number of sequence fragments
			_polap_log3_pipe "gfatools view \
        -S ${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa \
		    2>${_polap_output_dest} |
        grep '^S' | wc -l >${_polap_var_oga_summary}/${_pread_sel}/${i}.fragments"

			# Calculate the total number of bases in the sequence fragments
			_polap_log3_pipe "gfatools view \
        -S ${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa \
		    2>${_polap_output_dest} |
        grep '^S' | grep -o 'LN:i:[0-9]*' | awk -F: '{sum += \$3} END {print sum}' \
        >${_polap_var_oga_summary}/${_pread_sel}/${i}.bases"

			# Calculate the average depth
			cat "${_polap_var_oga_summary}/${_pread_sel}/${i}.gfa" |
				awk '/dp:i:/ {match($0, /dp:i:([0-9]+)/, arr); if (arr[1] != "") {sum += arr[1]; count++}} END {if (count > 0) print sum / count}' \
					>"${_polap_var_oga_summary}/${_pread_sel}/${i}.depth"

		else
			_polap_log1 "No such graph: ${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa"
		fi
	done

	# _polap_log1 NEXT: $0 flye2 -o "${_arg_outdir}" -j "${_arg_jnum}" -t "${_arg_threads}"
	if [[ "${_arg_menu[0]}" == "test-reads" ]]; then
		_polap_log1 "NEXT: $0 best-reads ..."
	else
		_polap_log1 "NEXT: $0 best-flye ..."
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
#
################################################################################
function _run_polap_select-reads { # selects reads mapped on a genome assembly
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
# Test organelle-assembly results with a chosen w value multiple times.
#
# Arguments:
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   -w ${_arg_single_min}: minimum minimap2 alignment length for a single contig
# Inputs:
#   ${_polap_var_mtcontigname}
#   ${_polap_var_ga_contigger_edges_fasta}
#   ${_polap_var_ga_contigger_edges_gfa} if no such file: ${_polap_var_ga_contigger_edges_fasta}
#   input long read data to use (priority order): 
#     1. ${_polap_var_outdir_lk_fq_gz}
#     2. ${_polap_var_outdir_nk_fq_gz}
#     3. ${_arg_long_reads}
# Outputs:
#   ${_polap_var_oga_contig}: map-reads output
#   ${_polap_var_oga_reads}: mapped read names
#   ${_polap_var_oga_seeds}: mapped reads
#   ${_polap_var_oga_subsample}: subsample of the mapped reads
#   ${_polap_var_oga_flye}: assemblies
#   ${_polap_var_oga_summary}: summary
#   ${_polap_var_oga_plot}: plot or table using range in contig folder and summary
# Menus:
#   ptgaul-reads [number of repeats]
#   intra-reads [number of repeats]
#   polap-reads [number of repeats]
Example: $(basename "$0") ${_arg_menu[0]} [ptgaul-reads] -w 3000
Example: $(basename "$0") ${_arg_menu[0]} intra-reads -w 5000
Example: $(basename "$0") ${_arg_menu[0]} polap-reads -w 3000 5
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	if [[ "${_arg_menu[1]}" == "view" ]]; then
		# use the summary to plot the number of bases
		# use the summary to plot the number of fragments
		_run_polap_test-reads

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	if [[ "${_arg_polap_reads}" == "on" ]]; then
		_arg_menu[1]="polap-reads"
		_polap_log0 "  use --polap-reads option, so menu1 becomes ${_arg_menu[1]}"
	fi

	case "${_arg_menu[1]}" in
	infile | ptgaul-reads)
		_arg_menu[1]="ptgaul-reads"
		_polap_log0 "  default set to read-selection: ${_arg_menu[1]}"
		;;
	polap-reads)
		_polap_log0 "  read-selection: ${_arg_menu[1]}"
		;;
	bridge-reads)
		_polap_log0 "  read-selection: ${_arg_menu[1]}"
		;;
	intra-reads)
		_polap_log0 "  read-selection: ${_arg_menu[1]}"
		;;
	*)
		_polap_log0 "ERROR: no such menu2 for select-reads: ${_arg_menu[1]}"
		return $RETURN_FAIL
		;;
	esac

	local _n_repeats="1"
	if [[ "${_arg_menu[2]}" != "outfile" ]]; then
		if [[ "${_n_repeats}" =~ ^[0-9]+$ ]]; then
			local _n_repeats="${_arg_menu[2]}"
		else
			_polap_log0 "ERROR: menu3 must be a number: ${_arg_menu[2]}"
			return 0
		fi
	fi
	_arg_select_read_range="${_arg_single_min},${_arg_single_min},${_n_repeats}"
	# _arg_flye="off"
	_run_polap_test-reads

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Best flye2
################################################################################
function _run_polap_flye2 { # executes Flye for an organelle-genome assembly
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	if [[ "${_arg_menu[1]}" == "infile" ]]; then
		_arg_menu[1]="ptgaul-reads"
		_polap_log0 "  default set to read-selection: ${_arg_menu[1]}"
	else
		_polap_log0 "  read-selection: ${_arg_menu[1]}"
	fi
	if [[ "${_arg_polap_reads}" == "on" ]]; then
		_arg_menu[1]="polap-reads"
		_polap_log0 "  use --polap-reads option, so menu1 becomes ${_arg_menu[1]}"
	fi
	local _pread_sel=${_arg_menu[1]}

	help_message=$(
		cat <<HEREDOC
# Execute Flye for an organelle-genome assembly.
#
# Arguments:
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   -t ${_arg_threads}: the number of CPU cores
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${_polap_var_oga_subsample}/${_pread_sel}/0.fq.gz
#   ${_polap_var_oga_contig}/contig.fa
# Outputs:
#   ${_polap_var_oga_assembly_graph_gfa}
#   ${_polap_var_oga_contigger_edges_gfa}
Example: $(basename "$0") ${_arg_menu[0]} [ptgaul-reads] -j <arg>
Example: $(basename "$0") ${_arg_menu[0]} intra-reads -j <arg>
Example: $(basename "$0") ${_arg_menu[0]} polap-reads -j <arg>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	if [[ "${_arg_redo}" == "on" ]]; then
		_polap_log3_cmd rm -rf "${_polap_var_oga}"/{00-assembly,10-consensus,20-repeat,30-contigger,40-polishing}
	fi

	_polap_log0 "flye2 runs on the assembly ${_arg_jnum}"
	local _contig_fa="${_polap_var_oga_contig}/contig.fa"
	local _long_reads="${_polap_var_oga_subsample}/${_pread_sel}/0.fq.gz"
	_polap_log1 "  input1: ${_polap_var_oga_contig}/contig.fa"
	_polap_log1 "  input2: ${_polap_var_oga_subsample}/${_pread_sel}/0.fq.gz"

	if [ ! -s "${_contig_fa}" ]; then
		_polap_log0 "ERROR: no selected-contig file: ${_contig_fa}"
		return
	fi

	if [[ ! -s "${_long_reads}" ]]; then
		_polap_log0 "ERROR: no long-read file: ${_long_reads}"
		return
	fi

	local _CONTIG_LENGTH=$(<"${_polap_var_oga_contig}/contig_total_length.txt")
	_polap_log2_cat "${_polap_var_oga_contig}/contig_total_length.txt"
	local _contig_length_bp=$(_polap_utility_convert_bp ${_CONTIG_LENGTH})

	_polap_log1 "  organelle genome size based on the seed contig selection: ${_contig_length_bp}"

	_polap_log1 "  executing the organelle-genome assembly using flye on ${_arg_jnum} ..."
	_polap_log1 "    input1: ${_long_reads}"
	_polap_log1 "    output1: ${_polap_var_oga}/30-contigger/graph_final.gfa"
	_polap_log1 "    output2: ${_polap_var_oga}/assembly_graph.gfa"
	# if [[ "${_arg_plastid}" == "on" ]]; then
	# 	_CONTIG_LENGTH=$((_CONTIG_LENGTH * 3))
	# fi
	local _command1="flye \
    ${_arg_flye_data_type} \
    ${_long_reads} \
		--out-dir ${_polap_var_oga} \
		--threads ${_arg_threads} \
		--asm-coverage ${_arg_flye_asm_coverage} \
		--genome-size ${_CONTIG_LENGTH} \
		2>${_polap_output_dest}"
	_polap_log3_pipe "${_command1}"

	rm -f "${_polap_var_output_oga_gfa}"
	ln -sf "${_polap_var_oga_assembly_graph_gfa}" "${_polap_var_output_oga_gfa}"
	_polap_log0 "  output: the assembly graph: ${_polap_var_oga_contigger_edges_gfa}"
	_polap_log0 "  output: the assembly graph: ${_polap_var_oga_assembly_graph_gfa}"
	_polap_log0 "  output: the assembly graph: ${_polap_var_output_oga_gfa}"
	jnum_next=$((_arg_jnum + 1))
	_polap_log1 "  create and edit ${_arg_outdir}/${_arg_jnum}/mt.contig.name-${jnum_next} and rerun assemble2"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Finishes a Flye organelle-genome assembly upto flye-polishing step
#
# Polishes an organelle-genome assembly using long-reads.
# Note: use the same options as flye2 menu.
# Arguments:
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   -t ${_arg_threads}: the number of CPU cores
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${MTDIR}/contig.fa
#   ${MTSEEDSDIR}/2.fq.gz
#   ${MTDIR}/30-contigger
# Outputs:
#   ${MTDIR}/assembly_graph.gfa
################################################################################
function _run_polap_flye-polishing { # finish a Flye organelle-genome assembly upto flye-polishing step
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	MTDIR="${_arg_outdir}"/${_arg_jnum}
	MTSEEDSDIR="${MTDIR}"/seeds

	help_message=$(
		cat <<HEREDOC
# Finishes the Flye organelle-genome assembly.
# Polishes an organelle-genome assembly using long-reads.
# Note: use the same options as flye2 menu.
# Arguments:
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   -t ${_arg_threads}: the number of CPU cores
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${MTDIR}/contig.fa
#   ${MTSEEDSDIR}/3.fq.gz
#   ${MTDIR}/30-contigger
# Outputs:
#   ${MTDIR}/assembly_graph.gfa
Example: $(basename "$0") ${_arg_menu[0]} -j <arg>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	if [ ! -s "${MTDIR}/contig.fa" ]; then
		echoall "ERROR: no selected-contig file [${MTDIR}/contig.fa]"
		echoerr "SUGGESTION: select-reads"
		exit $EXIT_SUCCESS
	fi

	if [ ! -s "${MTSEEDSDIR}/3.fq.gz" ]; then
		echoall "ERROR: no long-read file [${MTSEEDSDIR}/3.fq.gz]"
		echoerr "SUGGESTION: select-reads"
		exit $EXIT_SUCCESS
	fi

	_arg_menu[1]="polishing"
	_run_polap_flye2

	# CONTIG_LENGTH=$(seqkit stats -Ta "${MTDIR}"/contig.fa | csvtk cut -t -f "sum_len" | csvtk del-header)
	# _polap_log1 "INFO: organelle genome size based on contig selection: $CONTIG_LENGTH"
	#
	# _polap_log1 "INFO: polishing the organelle-genome assembly using flye ... be patient!"
	# _polap_log0 "please, wait for Flye long-read polishing of the organelle-genome assembly on ${_arg_jnum} ..."
	# _polap_log3_pipe "flye --nano-raw ${MTSEEDSDIR}/3.fq.gz \
	# 	--out-dir ${MTDIR} \
	# 	--threads ${_arg_threads} \
	# 	--asm-coverage ${_arg_flye_asm_coverage} \
	# 	--genome-size $CONTIG_LENGTH \
	# 	--resume \
	# 	2>$_polap_output_dest"
	#
	# _polap_log0 "  output: the long-read polished assembly graph: $PWD/${MTDIR}/assembly_graph.gfa"
	# _polap_log0 "  DO: extract a draft organelle genome sequence (mt.0.fasta) from the polished assembly graph"

	# echo "column -t ${_arg_outdir}/assembly_info_organelle_annotation_count.txt"
	# echoall NEXT: $(basename $0) check-coverage [-p ${_arg_unpolished_fasta}]
	# _polap_log1 NEXT: "$(basename "$0")" prepare-polishing -a "${_arg_short_read1}" -b "${_arg_short_read2}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Reports the organelle-genome assembly results.
################################################################################
function _run_polap_report-assembly { # report an organelle-genome assembly result
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Reports the organelle-genome assembly results.
#
# Arguments:
#   -o ${_arg_outdir}: output folder for BioProject
# Inputs:
#   ${_arg_outdir}: output folder for BioProject
# Outputs:
Example: $(basename $0) ${_arg_menu[0]} [-o ${_arg_outdir}]
Example: report-assembly -o PRJDB10540a 2>&1 | tr '\n' '\t' | sed 's/\t$/\n/'
HEREDOC
	)

	# Set variables for file paths
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_log0 "${help_message}" && return

	_polap_log1_log "reporting the organelle-genome assembly at ${_arg_outdir} ..."

	if [ -d "${_arg_outdir}" ]; then
		_polap_log2_file "the main output folder: ${_arg_outdir}"
	else
		_polap_log2 "ERROR: no such output folder; use -o option"
		exit $EXIT_SUCCESS
	fi

	_polap_log0 $(cut -f1 "${_polap_var_project_txt}")
	_polap_log0 $(cut -f1 "${_polap_var_project_sra_long_read}")
	_polap_log0 $(cut -f1 "${_polap_var_project_sra_short_read}")
	_polap_log0 $(cut -f1 "${_polap_var_project_species}")
	_polap_log0 $(cut -f1 "${_polap_var_project_mtdna_fasta2_accession}")

	# wc -l "${_polap_var_wga}/mt.contig.name-"? | awk '$2 != "total" {print $1}' | head -5 >&2

	# for i in "${_arg_select_contig_numbers[@]}"; do
	# 	# Call the function corresponding to the current number (index is i-1)
	# 	_arg_inum="${i}"
	# 	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
	# 	_polap_log0 $(cat "${_polap_var_wga}/mt.contig.name-${i}" | wc -l)
	# 	_polap_log0_cat "${_polap_var_mtdna_compare}"
	# done

	# Array to store the names of the original files
	files=($(ls "${_polap_var_wga}/mt.contig.name-"?))

	# Temporary array to store the paths of unique files
	unique_files=()

	# Function to check if a file is unique
	is_unique() {
		for unique_file in "${unique_files[@]}"; do
			if cmp -s "$1" "$unique_file"; then
				echo "$unique_file"
				return 1 # Not unique
			fi
		done
		return 0 # Unique
	}

	_polap_log1 "Checking for unique files and their matches:"

	# Iterate over the files to find unique ones
	for i in "${_arg_select_contig_numbers[@]}"; do
		# Call the function corresponding to the current number (index is i-1)
		local FDIR="${_polap_var_wga}"
		_arg_jnum="${i}"
		file="$FDIR"/mt.contig.name-${_arg_jnum}

		unique_file=$(is_unique "$file")
		if [ $? -eq 0 ]; then
			# If unique, add it to the unique_files array
			unique_files+=("$file")
			echo "$file is unique."

			_polap_var_mtcontigname="$file"
			_arg_inum="${i}"
		else
			_polap_log1 "$file is the same as $unique_file."
			_arg_inum="${unique_file##*-}"
		fi
		source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
		_polap_log0 $(cat "${_polap_var_wga}/mt.contig.name-${i}" | wc -l)
		_polap_log0_cat "${_polap_var_mtdna_compare}"
	done

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

# deleted these or move to run-polap-function-oga-v0.2.6.sh
# function _run_polap_x-v0.3.7-collect-reads() { # replaced by select-reads
# function _run_polap_x-v0.2.6-select-reads() { # selects reads mapped on a genome assembly
# function _run_polap_x-select-reads { # selects reads mapped on a genome assembly
# function _run_polap_x-v0.2.6-flye2() { # executes Flye for an organelle-genome assembly