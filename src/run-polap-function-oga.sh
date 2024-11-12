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
# Ensure that the current script is sourced only once
source "$script_dir/run-polap-function-include.sh"
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

source "$script_dir/run-polap-function-utilities.sh"

function _polap_oga_determine-long-read-file() {
	local -n result_ref=$1

	if [[ "${_arg_long_reads_is}" == "off" ]]; then
		if [[ -s "${_polap_var_base_lk_fq_gz}" ]]; then
			_polap_log2 "    we utilize all available size-limited long-read data for our analysis: ${_polap_var_base_lk_fq_gz}"
			result_ref="${_polap_var_base_lk_fq_gz}"
		elif [[ -s "${_polap_var_base_nk_fq_gz}" ]]; then
			_polap_log2 "    we utilize the sampled and size-limited long-read data for our analysis: ${_polap_var_base_nk_fq_gz}"
			result_ref="${_polap_var_base_nk_fq_gz}"
		else
			die "ERROR: no such file: ${_polap_var_base_lk_fq_gz}, ${_polap_var_base_nk_fq_gz}"
		fi
	else
		_polap_log2 "    we utilize the long-read data supplied through the command-line option -l."
		result_ref="${_arg_long_reads}"
	fi
}

################################################################################
# create a dummy assembly for preparing seed contigs
################################################################################
function _run_polap_prepare-seeds() { # prepare seed contigs in a not usual way
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
# Prepare contig seeds from multiple assemblies or sequences
#
# Arguments:
#   -i number1,number2,...
#   -j $JNUM: destination Flye organelle assembly number
#
#   -f contig sequence file
# Inputs:
#   mulitple assemblies
#   sequence file
# Outputs:
#   a dummy assembly
Example: $0 ${_arg_menu[0]} -i 4,6 -j 7
Example: $0 ${_arg_menu[0]} -f seed_contigs.fa
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	_polap_log0 "preparing seed contigs ..."
	_polap_log1 "  assembly: $INUM (source) -> $JNUM (target) ..."
	# if [[ -d "${ODIR}/${JNUM}" ]]; then
	# 	_polap_log0 "ERROR: target assembly number ${JNUM} already exists."
	# 	return $RETURN_SUCCESS
	# fi
	mkdir -p "${ODIR}/${JNUM}"
	if [[ "${_arg_final_assembly}" == "mt.1.fa" ]]; then
		_polap_log0 "  -i must have two or more assembly numbers"

		# Initialize a counter
		local counter=1

		# Output file for concatenated FASTA
		local output_file="${ODIR}/${JNUM}/seed_contigs.fa"

		# Clear or create the output file
		>"$output_file"
		# Set IFS to a comma to split the string by commas
		IFS=','

		# Use a for loop to iterate over each number
		for number in ${_arg_inum}; do
			local fasta_file="${ODIR}/${number}/assembly.fasta"

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
		# cp "${_arg_final_assembly}" "${ODIR}/${JNUM}/seed_contigs.fa"

		# Input FASTA file
		input_file="${_arg_final_assembly}"
		# Output FASTA file
		output_file="${ODIR}/${JNUM}/seed_contigs.fa"

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
	local KNUM=$((JNUM + 1))
	seqkit seq -n -i "${output_file}" >"${_polap_var_oga}/mt.contig.name-${KNUM}"
	_polap_log0_cat "${_polap_var_oga}/mt.contig.name-${KNUM}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
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
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   --
#   -r ${_arg_pair_min}: minimum minimap2 alignment length for a pair of contigs
#   -x ${_arg_bridge_min}: minimum long-read length for connecting the pair of contigs
#   -w ${_arg_single_min}: minimum minimap2 alignment length for a single contig
#
# Inputs:
#   _polap_var_mtcontigname="${_polap_var_ga}/mt.contig.name-${JNUM}"
#   ${_polap_var_mtcontigname} <- $INUM and $JNUM
#   ${MTCONTIGNAME} <- $INUM and $JNUM
#   ${_polap_var_contigger_edges_fasta} <- $INUM
#   ${_polap_var_base_nk_fq_gz}
#   ${_polap_var_base_lk_fq_gz}
#   ${_polap_var_base_l_fq_gz}
#
# Outputs:
#   _polap_var_oga_reads="${_polap_var_oga}/02-reads"
#   ${_polap_var_oga_reads}/contig.fa
#   ${_polap_var_oga_reads}/contig.paf
#   ${_polap_var_oga_reads}/contig.tab
################################################################################
function _run_polap_map-reads() { # selects reads mapped on a genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
# Map long reads on a Flye genome assembly.
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number (or 0)
#   -j $JNUM: destination Flye organelle assembly number
#   -l ${_arg_long_reads}: long-read data default:${_polap_var_base_lk_fq_gz}
# Inputs:
#   ${_polap_var_mtcontigname}
#   ${_polap_var_contigger_edges_fasta}
#   ${_polap_var_base_lk_fq_gz} or ${_polap_var_base_nk_fq_gz}
# Outputs:
#   ${_polap_var_oga_contig}/contig.fa
#   ${_polap_var_oga_contig}/contig.paf
#   ${_polap_var_oga_contig}/contig.tab
Example: $0 ${_arg_menu[0]} -i ${INUM} -j ${JNUM}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	_polap_log0 "mapping long-read data on the seed contigs ..."
	_polap_log1 "  assembly: $INUM (source) -> $JNUM (target) ..."
	_polap_log1 "  input1: ${_polap_var_mtcontigname}"
	_polap_log1 "  input2: ${_polap_var_contigger_edges_fasta}"

	if [[ -s "${_polap_var_oga_contig}/contig.tab" ]] && [[ "${_arg_redo}" == "on" ]]; then
		_polap_log0 "  found: ${_polap_var_oga_reads}/contig.tab, so skipping mapping long-read data ..."
		return
	fi

	if [ ! -s "${_polap_var_mtcontigname}" ]; then
		_polap_log0 "ERROR: no such mt.contig.name file: ${_polap_var_mtcontigname}"
		exit $EXIT_ERROR
	fi

	if [ ! -s "${_polap_var_contigger_edges_fasta}" ]; then
		_polap_log0 "ERROR: no assembly fasta file: ${_polap_var_contigger_edges_fasta}"
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

	# _polap_log3_cmd ln -s "$PWD/${_polap_var_base_nk_fq_gz}" -t "${_polap_var_oga}"
	# _polap_log3_cmd ln -s "$PWD/${_polap_var_base_lk_fq_gz}" -t "${_polap_var_oga}"
	_polap_log1 "  extracts contig sequeces from the assembly: ${_polap_var_contigger_edges_fasta}"
	_polap_log2 "    input1: ${_polap_var_contigger_edges_fasta}"
	_polap_log2 "    input2: ${_polap_var_mtcontigname}"
	_polap_log2 "    output: ${_polap_var_oga_contig}/contig.fa"
	_polap_log3_pipe "seqkit grep \
    -f ${_polap_var_mtcontigname} \
		${_polap_var_contigger_edges_fasta} \
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
			_polap_log3_cmd bash "$script_dir/run-polap-sh-half-cut.sh" "${_polap_var_oga_contig}" "${_polap_var_oga_contig}/contig.fa" "${_polap_var_mtcontigname}"

			# cp "${_polap_var_oga}/contig.fa" "${_polap_var_oga}/contig.fa-backup"
			# seqkit fx2tab --length --name "${_polap_var_oga}"/contig.fa -o "${_polap_var_oga}"/contig.fa.len >/dev/null 2>&1
			# A=$(cut -f2 "${_polap_var_oga}"/contig.fa.len)
			# B=$(echo "scale=0; $A/2" | bc)
			# C=$((B + 1))
			# seqkit subseq -r 1:"$B" "${_polap_var_oga}"/contig.fa -o "${_polap_var_oga}"/c1.fa >/dev/null 2>&1
			# seqkit subseq -r "$C":"$A" "${_polap_var_oga}"/contig.fa -o "${_polap_var_oga}"/c2.fa >/dev/null 2>&1
			# cat "${_polap_var_oga}"/c?.fa | seqkit replace -p '.+' -r 'edge_{nr}' -o "${_polap_var_oga}"/contig.fa >/dev/null 2>&1
			# cp "${_polap_var_mtcontigname}" "${MTCONTIGNAME}"-backup
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
		_polap_log3_pipe "minimap2 -cx map-ont \
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
	_polap_log3_cmd bash "$script_dir/run-polap-sh-minimap2-paf2tab.sh" "${_arg_min_read_length}" "${_polap_var_oga_contig}/contig.paf" "${_polap_var_oga_contig}/contig.tab"

	_polap_log1 "NEXT: $0 reads -o ${ODIR} -i ${INUM} -j ${JNUM}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
	return

	# put the backup to the original
	if [[ ${_arg_circularize} == "on" ]]; then
		_polap_log1 "  circularize option: putting seed contigs back to the original mt.contig.name and ${_polap_var_oga}/contig.fa"
		if [[ -s "${_polap_var_mtcontigname}"-backup ]]; then
			mv "${_polap_var_mtcontigname}"-backup "${MTCONTIGNAME}"
			mv "${_polap_var_oga}/contig.fa-backup" "${_polap_var_oga}/contig.fa"
		else
			echo "DEV: not implemented yet"
			exit $EXIT_ERROR
		fi
	fi

	_polap_log1 NEXT: $0 flye2 -o "$ODIR" -j "$JNUM" -t "${_arg_threads}" -c "${_arg_coverage}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Test organelle assemblies on a range of w values.
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   -s, --select-read-range <start,end,count>
#   --start-index <index>: to start at somewhere not start
# Inputs:
#   ${_polap_var_mtcontigname}
#   ${_polap_var_contigger_edges_fasta}
#   input long read data to use (priority order):
#     1. ${_polap_var_base_lk_fq_gz}
#     2. ${_polap_var_base_nk_fq_gz}
#     3. ${_arg_long_reads}
# Outputs:
#   ${_polap_var_oga_contig}: map-reads output
#   ${_polap_var_oga_reads}: mapped read names
#   ${_polap_var_oga_seeds}: mapped reads
#   ${_polap_var_oga_sample}: subsample of the mapped reads
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
function _run_polap_test-reads() { # selects reads mapped on a genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
# Test organelle assemblies on a range of w values.
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   -s, --select-read-range <start,end,count>
#   --start-index <index>: to start at somewhere not start
# Inputs:
#   ${_polap_var_mtcontigname}
#   ${_polap_var_contigger_edges_fasta}
#   input long read data to use (priority order): 
#     1. ${_polap_var_base_lk_fq_gz}
#     2. ${_polap_var_base_nk_fq_gz}
#     3. ${_arg_long_reads}
# Outputs:
#   ${_polap_var_oga_contig}: map-reads output
#   ${_polap_var_oga_reads}: mapped read names
#   ${_polap_var_oga_seeds}: mapped reads
#   ${_polap_var_oga_sample}: subsample of the mapped reads
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
Example: $0 ${_arg_menu[0]} [ptgaul-intra-base-length] --select-read-range 3000,27000,5
Example: $0 ${_arg_menu[0]} single-intra-base-length -i 1 -j 2
Example: $0 ${_arg_menu[0]} polap-rw-base-length --select-read-range 3000,27000,5
Example: $0 ${_arg_menu[0]} bridge-inter-base-length
Example: $0 ${_arg_menu[0]} polap-rw-base-length --select-read-range 3000,27000,5 --start-index 3
Example: $0 ${_arg_menu[0]} view polap-rw-base-length -i 2 -j 3
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	if [[ "${_arg_menu[1]}" == "preview" ]]; then
		seqkit stats -Ta "${_polap_var_oga_contig}/contig.fa" >&3

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
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

		cat ${_polap_var_oga_summary}/${_pread_sel}/*.bases \
			>"${_polap_var_oga_summary}/${_pread_sel}/bases.txt"
		cat ${_polap_var_oga_summary}/${_pread_sel}/*.fragments \
			>"${_polap_var_oga_summary}/${_pread_sel}/fragments.txt"
		cat ${_polap_var_oga_summary}/${_pread_sel}/*.depth \
			>"${_polap_var_oga_summary}/${_pread_sel}/depth.txt"

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
		_polap_log3_pipe "Rscript ${script_dir}/run-polap-r-test-reads-bar-graph.R \
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
		[ "$DEBUG" -eq 1 ] && set +x
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
			_polap_log3_pipe "Rscript ${script_dir}/run-polap-r-test-reads-bar-graph.R \
			-i ${output_file} \
			--out ${_polap_var_oga_plot}/${_pread_sel}-summary.pdf \
			>${_polap_output_dest} 2>&1"
		else
			_polap_log3_pipe "Rscript ${script_dir}/run-polap-r-test-reads-bar-graph.R \
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
		[ "$DEBUG" -eq 1 ] && set +x
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
	if [[ -d "${_polap_var_oga_reads}/${_pread_sel}" ]] && [[ "${_arg_start_index}" -eq 0 ]]; then
		if confirm "Do you want to redo the test-reads of ${_pread_sel} on assembly ${JNUM}"?; then
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
		rm -rf "${_polap_var_oga_sample}/${_pread_sel}"
		rm -rf "${_polap_var_oga_flye}/${_pread_sel}"
		rm -rf "${_polap_var_oga_summary}/${_pread_sel}"
		rm -rf "${_polap_var_oga_plot}/${_pread_sel}"
	fi

	_polap_log2 "  creating folders for read selection type: ${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_reads}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_seeds}/${_pread_sel}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_sample}/${_pread_sel}"
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
			_polap_log3_pipe "Rscript --vanilla ${script_dir}/run-polap-r-pairs.R \
	    -m ${_polap_var_mtcontigname} \
		  -t ${_polap_var_oga_contig}/contig.tab \
		  --out ${_polap_var_oga_reads}/${_pread_sel}/${i} \
      -w ${_arg_single_min} \
		  -r ${_arg_pair_min} \
		  -x ${_arg_bridge_min} \
      --all \
		  >${_polap_output_dest} 2>&1"
		else
			_polap_log3_pipe "Rscript --vanilla ${script_dir}/run-polap-r-bridge.R \
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

		if [[ -d "${_polap_var_base_lk_fq_gz}.split" ]]; then
			mkdir -p "${_polap_var_oga_seeds}/${_pread_sel}/${i}"
			ls "${_polap_var_base_lk_fq_gz}.split"/*.fq.gz |
				parallel seqtk subseq {} "${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.names" ">" "${_polap_var_oga_seeds}/${_pread_sel}/${i}/${_read_names}".{/}.fq

			# cat "${_polap_var_base_lk_fq_gz}.split"/*.fq.gz.fq |
			cat "${_polap_var_oga_seeds}/${_pread_sel}/${i}/${_read_names}".*.fq.gz.fq |
				gzip >"${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"
			rm -rf "${_polap_var_oga_seeds}/${_pread_sel}/${i}"
			# rm -f "${_polap_var_oga_seeds}/${_pread_sel}/${i}/${_read_names}".*.fq.gz.fq
			# rm -f "${_polap_var_base_lk_fq_gz}.split"/*.fq.gz.fq
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
		if [[ "$_expected_organelle_coverage" -gt "${_arg_coverage}" ]]; then
			if [[ "${_arg_coverage_check}" == "on" ]]; then
				local _rate=$(echo "scale=5; ${_arg_coverage}/$_expected_organelle_coverage" | bc)
				_polap_log0 "  long-read data reduction by rate of ${_rate} <= COV[${_arg_coverage}] / long-read organelle coverage[$_expected_organelle_coverage]"
				_polap_log1 "    sampling long-read data by ${_rate} ... wait ..."
				local _random_seed=${_arg_random_seed:-$RANDOM}
				# local _random_seed=11
				_polap_log1 "    random seed for reducing long reads mapped on potential seed contigs: ${_random_seed}"
				_polap_log3_pipe "seqkit sample \
        -p ${_rate} \
        -s ${_random_seed} \
			  ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz \
        -o ${_polap_var_oga_sample}/${_pread_sel}/${i}.fq.gz \
        2>${_polap_output_dest}"
				touch ${_polap_var_oga_sample}/${_pread_sel}/${i}.random.seed.${_random_seed}
			else
				_polap_log0 "    no reduction of the long-read data because of the option --no-coverage-check: expected coverage: ${_expected_organelle_coverage}"
				_polap_log3_cmd ln -s $(realpath "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz") "${_polap_var_oga_sample}/${_pread_sel}/${i}.fq.gz"
			fi
		else
			_polap_log0 "    no reduction of the long-read data because $_expected_organelle_coverage < ${_arg_coverage}"
			_polap_log3_cmd ln -s $(realpath "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz") "${_polap_var_oga_sample}/${_pread_sel}/${i}.fq.gz"
		fi

		_polap_log1 "  flye assembly for ${_pread_sel}"
		_polap_log2 "    input1: ${_polap_var_oga_sample}/${_pread_sel}/${i}.fq.gz"
		_polap_log2 "    output: ${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa"
		local _command1="flye \
      --nano-raw ${_polap_var_oga_sample}/${_pread_sel}/${i}.fq.gz \
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

	# _polap_log1 NEXT: $0 flye2 -o "$ODIR" -j "$JNUM" -t "${_arg_threads}" -c "${_arg_coverage}"
	if [[ "${_arg_menu[0]}" == "test-reads" ]]; then
		_polap_log1 "NEXT: $0 best-reads ..."
	else
		_polap_log1 "NEXT: $0 best-flye ..."
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
#
################################################################################
function _run_polap_select-reads() { # selects reads mapped on a genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
# Test organelle-assembly results with a chosen w value multiple times.
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   -w ${_arg_single_min}: minimum minimap2 alignment length for a single contig
# Inputs:
#   ${_polap_var_mtcontigname}
#   ${_polap_var_contigger_edges_fasta}
#   input long read data to use (priority order): 
#     1. ${_polap_var_base_lk_fq_gz}
#     2. ${_polap_var_base_nk_fq_gz}
#     3. ${_arg_long_reads}
# Outputs:
#   ${_polap_var_oga_contig}: map-reads output
#   ${_polap_var_oga_reads}: mapped read names
#   ${_polap_var_oga_seeds}: mapped reads
#   ${_polap_var_oga_sample}: subsample of the mapped reads
#   ${_polap_var_oga_flye}: assemblies
#   ${_polap_var_oga_summary}: summary
#   ${_polap_var_oga_plot}: plot or table using range in contig folder and summary
# Menus:
#   ptgaul-reads [number of repeats]
#   intra-reads [number of repeats]
#   polap-reads [number of repeats]
Example: $0 ${_arg_menu[0]} [ptgaul-reads] -w 3000
Example: $0 ${_arg_menu[0]} intra-reads -w 5000
Example: $0 ${_arg_menu[0]} polap-reads -w 3000 5
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
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		return
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
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Best flye2
################################################################################
function _run_polap_flye2() { # executes Flye for an organelle-genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh"

	if [[ "${_arg_menu[1]}" == "infile" ]]; then
		_arg_menu[1]="ptgaul-reads"
		_polap_log0 "  default set to read-selection: ${_arg_menu[1]}"
	else
		_polap_log0 "  read-selection: ${_arg_menu[1]}"
	fi
	local _pread_sel=${_arg_menu[1]}

	help_message=$(
		cat <<HEREDOC
# Execute Flye for an organelle-genome assembly.
#
# Arguments:
#   -j $JNUM: destination Flye organelle assembly number
#   -t ${_arg_threads}: the number of CPU cores
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${_polap_var_oga_sample}/${_pread_sel}/0.fq.gz
#   ${_polap_var_oga_contig}/contig.fa
# Outputs:
#   ${_polap_var_oga_assembly_graph_gfa}
#   ${_polap_var_oga_contigger_edges_gfa}
Example: $0 ${_arg_menu[0]} [ptgaul-reads] -j <arg>
Example: $0 ${_arg_menu[0]} intra-reads -j <arg>
Example: $0 ${_arg_menu[0]} polap-reads -j <arg>
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	if [[ "${_arg_redo}" == "on" ]]; then
		_polap_log3_cmd rm -rf "${_polap_var_oga}"/{00-assembly,10-consensus,20-repeat,30-contigger,40-polishing}
	fi

	_polap_log0 "flye2 runs on the assembly ${JNUM}"
	local _contig_fa="${_polap_var_oga_contig}/contig.fa"
	local _long_reads="${_polap_var_oga_sample}/${_pread_sel}/0.fq.gz"
	_polap_log1 "  input1: ${_polap_var_oga_contig}/contig.fa"
	_polap_log1 "  input2: ${_polap_var_oga_sample}/${_pread_sel}/0.fq.gz"

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

	_polap_log1 "  executing the organelle-genome assembly using flye on $JNUM ..."
	_polap_log1 "    input1: ${_long_reads}"
	_polap_log1 "    output1: ${_polap_var_oga}/30-contigger/graph_final.gfa"
	_polap_log1 "    output2: ${_polap_var_oga}/assembly_graph.gfa"
	local _command1="flye \
    --nano-raw ${_long_reads} \
		--out-dir ${_polap_var_oga} \
		--threads ${_arg_threads} \
		--asm-coverage ${_arg_flye_asm_coverage} \
		--genome-size ${_CONTIG_LENGTH} \
		2>${_polap_output_dest}"
	_polap_log3_pipe "${_command1}"

	_polap_log0 "  output: the assembly graph: ${_polap_var_oga_contigger_edges_gfa}"
	_polap_log0 "  output: the assembly graph: ${_polap_var_oga_assembly_graph_gfa}"
	jnum_next=$((JNUM + 1))
	_polap_log1 "  create and edit $ODIR/$JNUM/mt.contig.name-${jnum_next} and rerun assemble2"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Finishes a Flye organelle-genome assembly upto flye-polishing step
#
# Polishes an organelle-genome assembly using long-reads.
# Note: use the same options as flye2 menu.
# Arguments:
#   -j $JNUM: destination Flye organelle assembly number
#   -t ${_arg_threads}: the number of CPU cores
#   -c ${_arg_coverage}: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${MTDIR}/contig.fa
#   ${MTSEEDSDIR}/2.fq.gz
#   ${MTDIR}/30-contigger
# Outputs:
#   ${MTDIR}/assembly_graph.gfa
################################################################################
function _run_polap_flye-polishing() { # finish a Flye organelle-genome assembly upto flye-polishing step
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	MTDIR="$ODIR"/$JNUM
	MTSEEDSDIR="${MTDIR}"/seeds

	help_message=$(
		cat <<HEREDOC
# Finishes the Flye organelle-genome assembly.
# Polishes an organelle-genome assembly using long-reads.
# Note: use the same options as flye2 menu.
# Arguments:
#   -j $JNUM: destination Flye organelle assembly number
#   -t ${_arg_threads}: the number of CPU cores
#   -c ${_arg_coverage}: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${MTDIR}/contig.fa
#   ${MTSEEDSDIR}/3.fq.gz
#   ${MTDIR}/30-contigger
# Outputs:
#   ${MTDIR}/assembly_graph.gfa
Example: $0 ${_arg_menu[0]} -j <arg>
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
	# _polap_log0 "please, wait for Flye long-read polishing of the organelle-genome assembly on $JNUM ..."
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

	# echo "column -t $ODIR/assembly_info_organelle_annotation_count.txt"
	# echoall NEXT: $(basename $0) check-coverage [-p ${_arg_unpolished_fasta}]
	# _polap_log1 NEXT: "$(basename "$0")" prepare-polishing -a "${_arg_short_read1}" -b "${_arg_short_read2}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Reports the organelle-genome assembly results.
################################################################################
function _run_polap_report-assembly() { # report an organelle-genome assembly result
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	# Help message
	local help_message=$(
		cat <<HEREDOC
# Reports the organelle-genome assembly results.
#
# Arguments:
#   -o ${ODIR}: output folder for BioProject
# Inputs:
#   ${ODIR}: output folder for BioProject
# Outputs:
Example: $(basename $0) ${_arg_menu[0]} [-o ${ODIR}]
Example: report-assembly -o PRJDB10540a 2>&1 | tr '\n' '\t' | sed 's/\t$/\n/'
HEREDOC
	)

	# Set variables for file paths
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_log0 "${help_message}" && return

	_polap_log1_log "reporting the organelle-genome assembly at ${ODIR} ..."

	if [ -d "${ODIR}" ]; then
		_polap_log2_file "the main output folder: ${ODIR}"
	else
		_polap_log2 "ERROR: no such output folder; use -o option"
		exit $EXIT_SUCCESS
	fi

	_polap_log0 $(cut -f1 "${_polap_var_bioproject_txt}")
	_polap_log0 $(cut -f1 "${_polap_var_bioproject_sra_long_read}")
	_polap_log0 $(cut -f1 "${_polap_var_bioproject_sra_short_read}")
	_polap_log0 $(cut -f1 "${_polap_var_bioproject_species}")
	_polap_log0 $(cut -f1 "${_polap_var_bioproject_mtdna_fasta2_accession}")

	# wc -l "${_polap_var_wga}/mt.contig.name-"? | awk '$2 != "total" {print $1}' | head -5 >&2

	# for i in "${_arg_select_contig_numbers[@]}"; do
	# 	# Call the function corresponding to the current number (index is i-1)
	# 	INUM="${i}"
	# 	source "$script_dir/polap-variables-common.sh" # '.' means 'source'
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
		JNUM="${i}"
		file="$FDIR"/mt.contig.name-$JNUM

		unique_file=$(is_unique "$file")
		if [ $? -eq 0 ]; then
			# If unique, add it to the unique_files array
			unique_files+=("$file")
			echo "$file is unique."

			MTCONTIGNAME="$file"
			INUM="${i}"
		else
			_polap_log1 "$file is the same as $unique_file."
			INUM="${unique_file##*-}"
		fi
		source "$script_dir/polap-variables-common.sh" # '.' means 'source'
		_polap_log0 $(cat "${_polap_var_wga}/mt.contig.name-${i}" | wc -l)
		_polap_log0_cat "${_polap_var_mtdna_compare}"
	done

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
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
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
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
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	# for contigs
	#	_polap_var_contigger_edges_fasta=o/30-contigger/contigs.fasta
	# for edges

	help_message=$(
		cat <<HEREDOC
# Select long reads mapped on a genome assembly.
#
# Arguments:
#   -i ${INUM}: source Flye (usually whole-genome) assembly number (or 0)
#   -j ${JNUM}: destination Flye organelle assembly number
#   -l ${_arg_long_reads}
# Inputs:
#   ${_polap_var_mtcontigname}
#   ${_polap_var_contigger_edges_fasta}
#   input long read data: 1. ${_polap_var_base_lk_fq_gz} <- input data used for the whole-genome assembly
#                         2. ${_arg_long_reads}          <- input long-read data
# Outputs:
#   ${_polap_var_oga_seeds}/ptgaul.fq.gz   <- ptGAUL read-selection
#   ${_polap_var_oga_seeds}/single.fq.gz   <- single-mapped reads, if we have many such reads
#   ${_polap_var_oga_seeds}/pair.fq.gz     <- only pair-mapped reads (not useful)
#   ${_polap_var_oga_seeds}/combined.fq.gz <- single-mapped reads plus pair-mapped reads
Example: $0 ${_arg_menu[0]} -i ${INUM} -j ${JNUM}
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
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	_polap_log0 "selecting long-reads mapped on the seed contigs in file"

	_polap_log1 "  please, wait for a long-read data selection ..."
	_polap_log1 "    $INUM (source) -> $JNUM (target) ..."
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
	_polap_log3_pipe "Rscript --vanilla ${script_dir}/run-polap-r-pairs.R \
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
		local _source_long_reads_fq="${_polap_var_base_lk_fq_gz}"
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

	_polap_log1 NEXT: $0 flye2 -o "$ODIR" -j "$JNUM" -t "${_arg_threads}" -c "${_arg_coverage}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Selects reads mapped on a genome assembly.
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   -r ${_arg_pair_min}: minimum minimap2 alignment length for a pair of contigs
#   -x ${_arg_bridge_min}: minimum long-read length for connecting the pair of contigs
#   -w ${_arg_single_min}: minimum minimap2 alignment length for a single contig
# Inputs:
#   ${MTCONTIGNAME}
#   ${_polap_var_contigger_edges_fasta}
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
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh" # '.' means 'source'
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	local MTDIR="${_polap_var_oga}"                              # target: ${_polap_var_oga}
	local MTSEEDSDIR="${_polap_var_oga}/seeds"                   # ${_polap_var_seeds} for oga-class
	local MTCONTIGNAME="${_polap_var_ga}/mt.contig.name-${JNUM}" # ${_polap_var_mtcontigname}

	# for contigs
	#	_polap_var_contigger_edges_fasta=o/30-contigger/contigs.fasta
	# for edges

	help_message=$(
		cat <<HEREDOC
# Select long reads mapped on a genome assembly.
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number (or 0)
#   -j $JNUM: destination Flye organelle assembly number
#   --rwx ${_arg_rwx}: set the option values of -r, -x, -w to the same one
# Inputs:
#   ${MTCONTIGNAME}
#   ${_polap_var_contigger_edges_fasta}
# Outputs:
#   ${MTDIR}/contig.fa
#   ${MTSEEDSDIR}/1.fq.gz <- single-mapped reads
#   ${MTSEEDSDIR}/2.fq.gz <- reduced single-mapped reads, if we have many such reads
#   ${MTSEEDSDIR}/3.fq.gz <- reduced single-mapped reads plus pair-mapped reads
Example: $0 ${_arg_menu[0]} -i ${INUM} -j ${JNUM} --rwx ${_arg_rwx}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	_polap_log0 "selecting long-reads mapped on the seed contigs in file: ${MTCONTIGNAME} and ${_polap_var_contigger_edges_fasta}"
	_polap_log1 "  input1: ${MTCONTIGNAME}"
	_polap_log1 "  input2: ${_polap_var_contigger_edges_fasta}"

	if [ ! -s "${MTCONTIGNAME}" ]; then
		_polap_log0 "ERROR: no such mt.contig.name file: ${MTCONTIGNAME}"
		exit $EXIT_ERROR
	fi

	if [ ! -s "${_polap_var_contigger_edges_fasta}" ]; then
		_polap_log0 "ERROR: no assembly fasta file: ${_polap_var_contigger_edges_fasta}"
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

	_polap_log1 "  creating ${ODIR}/organelle-assembly_${INUM}-${JNUM}"
	echo "$CMD" >"${ODIR}/organelle-assembly_${INUM}-${JNUM}"

	_polap_log1 "  creates ${MTSEEDSDIR}"
	_polap_log3_cmd mkdir -p "${MTSEEDSDIR}"
	_polap_log3_cmd ln -s "$PWD/${_polap_var_base_nk_fq_gz}" -t "${MTDIR}"
	_polap_log1 "  extracts contig sequeces from the assembly: ${_polap_var_contigger_edges_fasta}"
	_polap_log2 "    input1: ${_polap_var_contigger_edges_fasta}"
	_polap_log2 "    input2: ${MTCONTIGNAME}"
	_polap_log2 "    output: ${MTDIR}/contig.fa"
	_polap_log3_pipe "seqkit grep \
    --threads ${_arg_threads} \
    -f ${MTCONTIGNAME} \
		${_polap_var_contigger_edges_fasta} \
		-o ${MTDIR}/contig.fa \
		2>${_polap_output_dest}"

	# we could circularize a single contig.
	local contig_count=$(wc -l <"${MTCONTIGNAME}")
	_polap_log1 "  number of seed contigs: ${contig_count}"
	if [[ ${_arg_circularize} == "on" ]]; then
		if [ "$contig_count" -eq 1 ]; then

			_polap_log1 "  circularizes the single seed contig ..."
			_polap_log2 "    input1: ${MTDIR}"
			_polap_log2 "    input2: ${MTDIR}/contig.fa"
			_polap_log2 "    input3: ${MTCONTIGNAME}"
			_polap_log2 "    output1: ${MTCONTIGNAME}-backup"
			_polap_log2 "    output2: (new) ${MTCONTIGNAME}"
			_polap_log2 "    output3: ${MTDIR}/contig.fa-backup"
			_polap_log2 "    output4: (new) ${MTDIR}/contig.fa"
			_polap_log3_cmd bash "$script_dir/run-polap-sh-half-cut.sh" "${MTDIR}" "${MTDIR}/contig.fa" "${MTCONTIGNAME}"

			# cp "${MTDIR}/contig.fa" "${MTDIR}/contig.fa-backup"
			# seqkit fx2tab --length --name "${MTDIR}"/contig.fa -o "${MTDIR}"/contig.fa.len >/dev/null 2>&1
			# A=$(cut -f2 "${MTDIR}"/contig.fa.len)
			# B=$(echo "scale=0; $A/2" | bc)
			# C=$((B + 1))
			# seqkit subseq -r 1:"$B" "${MTDIR}"/contig.fa -o "${MTDIR}"/c1.fa >/dev/null 2>&1
			# seqkit subseq -r "$C":"$A" "${MTDIR}"/contig.fa -o "${MTDIR}"/c2.fa >/dev/null 2>&1
			# cat "${MTDIR}"/c?.fa | seqkit replace -p '.+' -r 'edge_{nr}' -o "${MTDIR}"/contig.fa >/dev/null 2>&1
			# cp "${MTCONTIGNAME}" "${MTCONTIGNAME}"-backup
			# echo -e "edge_1\nedge_2" >"${MTCONTIGNAME}"

			_polap_log1 "    creating new ${MTDIR}/contig.fa and ${MTCONTIGNAME}"

		else
			_polap_log0 "Not implemented yet!"
			exit $EXIT_ERROR
			# "$script_dir"/run-polap-single.R "${MTSEEDSDIR}"/contig.tab "${MTSEEDSDIR}" "${_arg_single_min}" >/dev/null 2>&1
			# cat "${MTSEEDSDIR}"/single.names | sort | uniq >"${MTSEEDSDIR}"/1.names
			# echo "INFO: creates long read single name in ${MTSEEDSDIR}/1.names"
		fi
	fi

	_polap_log1 "  please, wait for a long-read data selection ..."
	_polap_log1 "    $INUM (source) -> $JNUM (target) ..."
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
	_polap_log2 "    input2: ${_polap_var_base_nk_fq_gz}"
	_polap_log2 "    output: ${MTDIR}/contig.paf"
	if [[ -s "${MTDIR}"/contig.paf ]] && [[ "${_arg_redo}" = "off" ]]; then
		_polap_log1 "  found: ${MTDIR}/contig.paf, skipping the minimap2 mapping step ..."
	else
		_polap_log3_pipe "minimap2 -cx map-ont \
      ${MTDIR}/contig.fa \
      ${_polap_var_base_nk_fq_gz} \
      -t ${_arg_threads} \
      -o ${MTDIR}/contig.paf \
      >${_polap_output_dest} 2>&1"
	fi

	_polap_log1 "  filtering out reads less than ${_arg_min_read_length} bp and converting PAF to TAB ..."
	_polap_log2 "    input1: ${MTDIR}/contig.paf"
	_polap_log2 "    output: ${MTDIR}/contig.tab"
	# cut -f1-11 "${MTDIR}"/contig.paf | awk -v minlength="${_arg_min_read_length}" '{if ($2>=minlength) {print}}' >"${MTDIR}"/contig.tab
	_polap_log3_cmd bash "$script_dir/run-polap-sh-minimap2-paf2tab.sh" "${_arg_min_read_length}" "${MTDIR}/contig.paf" "${MTDIR}/contig.tab"

	# MT: MPAIR=3000 MBRIDGE=3000 MSINGLE=3000
	# PT: MPAIR=1000 MBRIDGE=5000 MSINGLE=0
	_polap_log1 "  selecting seeds that meet the selection criteria ..."
	_polap_log2 "    bridge (POLAP pairs bridge minimum)=${_arg_bridge_min}"
	_polap_log2 "    p_mapping (POLAP pairs alignment minimum)=${_arg_pair_min}"
	_polap_log2 "    s_mapping (POLAP single alignment minimum)=${_arg_single_min}"
	_polap_log2 "    min_len_read=${_arg_min_read_length}"
	_polap_log2 "    input1: ${MTCONTIGNAME}"
	_polap_log2 "    input2: ${MTDIR}/contig.tab"
	_polap_log2 "    output: ${MTSEEDSDIR}/single.names for reads mapped on a single contig"
	_polap_log2 "    output: ${MTSEEDSDIR}/<edge1>-<edge2>.name for reads mapped on contig pairs"
	_polap_log3_pipe "Rscript --vanilla ${script_dir}/run-polap-pairs.R \
		${MTCONTIGNAME} \
		${MTDIR}/contig.tab \
		${MTSEEDSDIR} \
		${_arg_pair_min} \
    ${_arg_bridge_min} \
    ${_arg_single_min} \
    >${_polap_output_dest} 2>&1"
	# "$script_dir"/run-polap-pairs.R "${MTCONTIGNAME}" ${MTDIR}/contig.tab ${MTSEEDSDIR} ${_arg_pair_min} ${_arg_bridge_min} ${_arg_single_min} >/dev/null 2>&1

	# cat "${MTSEEDSDIR}/"*".name" "${MTSEEDSDIR}"/single.names | sort | uniq >"${MTSEEDSDIR}"/1.names
	_polap_log1 "  creates names of reads mapped on a single name in ${MTSEEDSDIR}/1.names"
	_polap_log2 "    input1: ${MTSEEDSDIR}/single.names"
	_polap_log2 "    output: ${MTSEEDSDIR}/1.names"
	_polap_log3_pipe "cat ${MTSEEDSDIR}/single.names | sort | uniq >${MTSEEDSDIR}/1.names"

	# seqkit grep --threads ${_arg_threads} -f "${MTSEEDSDIR}"/1.names ${_polap_var_base_nk_fq_gz} -o "${MTSEEDSDIR}"/1.fq.gz >/dev/null 2>&1
	_polap_log1 "  creating reads mapped on a single contig in ${MTSEEDSDIR}/1.fq.gz"
	_polap_log2 "    input1: ${_polap_var_base_nk_fq_gz}"
	_polap_log2 "    input2: ${MTSEEDSDIR}/1.names"
	_polap_log2 "    output: ${MTSEEDSDIR}/1.fq.gz"
	seqtk subseq "${_polap_var_base_nk_fq_gz}" "${MTSEEDSDIR}"/1.names | gzip >"${MTSEEDSDIR}"/1.fq.gz

	_polap_log1 "  computing the total size of the single-mapped reads ..."
	_polap_log2 "    input1: ${MTSEEDSDIR}/1.fq.gz"
	local _contig_length_bp=$(_polap_utility_convert_bp ${CONTIG_LENGTH})
	_polap_log2 "    input2 (seed contig size): ${_contig_length_bp}"
	local TOTAL_LENGTH=$(seqkit stats -Ta "${MTSEEDSDIR}"/1.fq.gz | csvtk cut -t -f "sum_len" | csvtk del-header)
	local _total_length_bp=$(_polap_utility_convert_bp ${TOTAL_LENGTH})
	_polap_log2 "    result1 (total size of reads mapped on a single contig): ${_total_length_bp}"
	local EXPECTED_ORGANELLE_COVERAGE=$((TOTAL_LENGTH / CONTIG_LENGTH))
	_polap_log2 "    result2 (expected organelle coverage): ${EXPECTED_ORGANELLE_COVERAGE}x"

	_polap_log1 "  reducing the single-mapped reads upto the coverage of ${_arg_coverage}x"
	_polap_log2 "    input1: ${MTSEEDSDIR}/1.fq.gz"
	_polap_log2 "    output: ${MTSEEDSDIR}/2.fq.gz"
	if [[ "${_arg_test}" == "on" ]]; then
		_polap_log0 "    OPTION: --test : No reduction of the test long-read data"
		_polap_log3_cmd ln -s $(realpath "${MTSEEDSDIR}"/1.fq.gz) "${MTSEEDSDIR}"/2.fq.gz
	elif [[ "${_arg_coverage_check}" == "off" ]]; then
		_polap_log0 "    OPTION: --no-coverage-check : No reduction of the long-read data"
		_polap_log3_cmd ln -s $(realpath "${MTSEEDSDIR}"/1.fq.gz) "${MTSEEDSDIR}"/2.fq.gz
	else
		if [ "$EXPECTED_ORGANELLE_COVERAGE" -lt "${_arg_coverage}" ]; then
			_polap_log1 "    no reduction of the long-read data because $EXPECTED_ORGANELLE_COVERAGE < ${_arg_coverage}"
			_polap_log3_cmd ln -s $(realpath "${MTSEEDSDIR}"/1.fq.gz) "${MTSEEDSDIR}"/2.fq.gz
		else
			_polap_log1 "    SUGGESTION: you might want to increase the minimum read lengths (use --rwx or -m) because you have enough long-read data."
			RATE=$(echo "scale=10; ${_arg_coverage}/$EXPECTED_ORGANELLE_COVERAGE" | bc)
			_polap_log0 "  long-read data reduction by rate of $RATE <= COV[${_arg_coverage}] / long-read organelle coverage[$EXPECTED_ORGANELLE_COVERAGE]"
			_polap_log1 "    sampling long-read data by $RATE ... wait ..."
			# seqkit sample -p "$RATE" "${MTSEEDSDIR}/1.fq.gz" -o "${MTSEEDSDIR}/2.fq.gz" >/dev/null 2>&1
			local seed=${_arg_random_seed:-$RANDOM}
			_polap_log1 "    random seed for reducing the single-mapped long-read data: ${seed}"
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
	# seqtk subseq ${_polap_var_base_nk_fq_gz} ${MTSEEDSDIR}/1.names.2 | gzip >${MTSEEDSDIR}/3.fq.gz
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
		# seqkit grep --threads ${_arg_threads} -f "${MTSEEDSDIR}"/1.names.2 ${_polap_var_base_nk_fq_gz} -o "${MTSEEDSDIR}"/2.fq.gz >/dev/null 2>&1

		_polap_log1 "  extracting single- and pair-mapped reads ..."
		_polap_log2 "    input1: ${_polap_var_base_nk_fq_gz}"
		_polap_log2 "    input1: ${MTSEEDSDIR}/1.names.2"
		_polap_log2 "    output: ${MTSEEDSDIR}/3.fq.gz"
		_polap_log3_pipe "seqtk subseq ${_polap_var_base_nk_fq_gz} ${MTSEEDSDIR}/1.names.2 | gzip >${MTSEEDSDIR}/3.fq.gz"
	else
		_polap_log1 "    no pair-mapped reads: 2.fq.gz -> 3.fq.gz"
		_polap_log3_cmd ln -s $(realpath "${MTSEEDSDIR}"/2.fq.gz) "${MTSEEDSDIR}"/3.fq.gz
	fi
	_polap_log1 "  organelle reads in ${MTSEEDSDIR}/3.fq.gz"

	# put the backup to the original
	if [[ ${_arg_circularize} == "on" ]]; then
		_polap_log1 "  circularize option: putting seed contigs back to the original mt.contig.name and ${MTDIR}/contig.fa"
		if [[ -s "${MTCONTIGNAME}"-backup ]]; then
			mv "${MTCONTIGNAME}"-backup "${MTCONTIGNAME}"
			mv "${MTDIR}/contig.fa-backup" "${MTDIR}/contig.fa"
		else
			echo "DEV: not implemented yet"
			exit $EXIT_ERROR
		fi
	fi

	_polap_log1 NEXT: $0 flye2 -o "$ODIR" -j "$JNUM" -t "${_arg_threads}" -c "${_arg_coverage}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# simply use read names to extract reads from a long-read data file
################################################################################
function _run_polap_x-select-reads() { # selects reads mapped on a genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	local MTDIR="${_polap_var_oga}"                              # target: ${_polap_var_oga}
	local MTSEEDSDIR="${_polap_var_oga}/seeds"                   # ${_polap_var_seeds} for oga-class
	local MTCONTIGNAME="${_polap_var_ga}/mt.contig.name-${JNUM}" # ${_polap_var_mtcontigname}

	# for contigs
	#	_polap_var_contigger_edges_fasta=o/30-contigger/contigs.fasta
	# for edges

	help_message=$(
		cat <<HEREDOC
# Extract long reads mapped on a genome assembly.
#
# Arguments:
#   -j ${JNUM}: destination Flye organelle assembly number
#   -l ${_arg_long_reads}
# Inputs:
#   input long read data: 1. ${_polap_var_base_nk_fq_gz} <- input data used for the whole-genome assembly
#                         2. ${_arg_long_reads}          <- input long-read data
#   ${ODIR}/${JNUM}/seeds/1.names         <- single-mapped read
#   ${ODIR}/${JNUM}/seeds/single.names.2  <- reduced single-mapped read
#   ${ODIR}/${JNUM}/seeds/1.names.2       <- reduced single-mapped read + pair-mapped read
#   ${_polap_var_contigger_edges_fasta}
# Outputs:
#   ${_polap_var_oga_seeds}/2.fq.gz <- reduced single-mapped reads, if we have many such reads
#   ${_polap_var_oga_seeds}/3.fq.gz <- reduced single-mapped reads plus pair-mapped reads
Example: $0 ${_arg_menu[0]} -j ${JNUM}
Example: $0 ${_arg_menu[0]} -j ${JNUM} -l l.fq
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
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	_polap_log0 "extracting reads using read names ..."

	if [[ "${_arg_long_reads_is}" == "off" ]]; then
		local _source_long_reads_fq="${_polap_var_base_nk_fq_gz}"
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

	_polap_log1 "  extracts contig sequeces from the assembly: ${_polap_var_contigger_edges_fasta}"
	_polap_log2 "    input1: ${_polap_var_contigger_edges_fasta}"
	_polap_log2 "    input2: ${MTCONTIGNAME}"
	_polap_log2 "    output: ${MTDIR}/contig.fa"
	_polap_log3_pipe "seqkit grep \
    --threads ${_arg_threads} \
    -f ${MTCONTIGNAME} \
		${_polap_var_contigger_edges_fasta} \
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
	# seqtk subseq ${_polap_var_base_nk_fq_gz} ${MTSEEDSDIR}/1.names.2 | gzip >${MTSEEDSDIR}/3.fq.gz

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
			_polap_log1 "  executing the organelle-genome long-read flye polishing assembly on $JNUM ..."
			_polap_log1 "    input1: ${_polap_var_oga_seeds}/${_i_seed}.fq.gz"
			_polap_log1 "    output: ${_polap_var_oga_seeds}/${_i_seed}/assembly_graph.gfa"
		else
			_polap_log1 "  executing the organelle-genome assembly using flye on $JNUM ..."
			_polap_log1 "    input1: ${_polap_var_oga_seeds}/${_i_seed}.fq.gz"
			_polap_log1 "    output: ${_polap_var_oga_seeds}/${_i_seed}/30-contigger/graph_final.gfa"
		fi

		local _command1="flye \
    --nano-raw ${_polap_var_oga_seeds}/${_i_seed}.fq.gz \
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

	_polap_log1 NEXT: $0 flye2 -o "$ODIR" -j "$JNUM" -t "${_arg_threads}" -c "${_arg_coverage}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Executes Flye for an organelle-genome assembly
#
# Arguments:
#   -j $JNUM: destination Flye organelle assembly number
#   -t ${_arg_threads}: the number of CPU cores
#   -c ${_arg_coverage}: the Flye's coverage option
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
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh"

	local MTDIR="$ODIR"/$JNUM
	MTSEEDSDIR="${MTDIR}"/seeds

	help_message=$(
		cat <<HEREDOC
# Executes Flye for an organelle-genome assembly
# Arguments:
#   -j $JNUM: destination Flye organelle assembly number
#   -t ${_arg_threads}: the number of CPU cores
#   -c ${_arg_coverage}: the Flye's coverage option
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
Example: $0 ${_arg_menu[0]} -j <arg>
Example: $0 ${_arg_menu[0]} -j <arg> polishing
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
		_polap_log1 "  executing the organelle-genome long-read flye polishing assembly on $JNUM ..."
		_polap_log1 "    input1: ${MTSEEDSDIR}/3.fq.gz"
		_polap_log1 "    output: ${MTDIR}/assembly_graph.gfa"
	else
		_polap_log1 "  executing the organelle-genome assembly using flye on $JNUM ..."
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
		_polap_log1 "NEXT: $0 prepare-polishing -a ${_arg_short_read1} -b ${_arg_short_read2}"
	else
		_polap_log0 "  output: the assembly fasta: ${_polap_var_oga_contigger_edges_fasta}"
		_polap_log0 "  output: the assembly graph: ${_polap_var_oga_contigger_edges_gfa}"
		jnum_next=$((JNUM + 1))
		_polap_log1 "  create and edit $ODIR/$JNUM/mt.contig.name-${jnum_next}"
		_polap_log1 "NEXT: $0 assemble2 -o ${ODIR} -i ${JNUM} -j ${jnum_next}"
		_polap_log1 or you could finish with Flye organelle-genome assembly with its polishing stage.
		_polap_log1 "NEXT: $0 flye-polishing -o ${ODIR} -j $${JNUM}"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
