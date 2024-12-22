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
#
################################################################################
function _run_polap_prepare-directional-seeds { #
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
This tool creates the seed contigs in FASTA format.

Inputs
------

- source assembly number: i
- destination assembly number: j
- i/mt.contig.name-j
- i/assembly_graph.gfa

Outputs
-------

- j/01-contig/contig1.fa
- j/01-contig/contig1.name.txt
- j/01-contig/contig2.fa
- j/01-contig/contig2.name.txt

Usage
-----
$(basename "$0") ${_arg_menu[0]} -i 1 -j 2
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	if [[ "${_arg_menu[1]}" == "preview" ]]; then
		_polap_log0_cat "${_polap_var_mtcontigname}"
		ls -l "${_polap_var_ga_assembly_graph_gfa}" >&3

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "view" ]]; then
		_polap_log0_cat "${_polap_var_oga_contig}/contig1.name.txt"
		seqkit stats "${_polap_var_oga_contig}/contig1.fa" >&3
		_polap_log0_cat "${_polap_var_oga_contig}/contig2.name.txt"
		seqkit stats "${_polap_var_oga_contig}/contig2.fa" >&3

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	_polap_log0 "preparing seed contigs ..."
	_polap_log1 "  assembly: ${_arg_inum} (source) -> ${_arg_jnum} (target) ..."

	if [[ "${_arg_redo}" == "on" ]]; then
		_polap_log3_cmd rm -rf "${_arg_outdir}/${_arg_jnum}"
	fi

	_polap_log3_cmd mkdir -p "${_polap_var_oga_contig}"

	# gfa -> fasta
	# i/30-contigger/graph_final.gfa
	# i/mt.contig.name-j
	#
	# j/01-contig/contig1.fa
	# j/01-contig/contig1.name.txt
	#
	# edge_number_sign -> contig_number
	# edge_1+, edge_2-, edge_7-, edge_2+ -> contig_1
	# edge_3+ -> contig_2
	#
	# j/01-contig/contig2.fa
	# j/01-contig/contig2.name.txt
	# edge_1+
	# edge_7-
	# edge_3+

	# Main 1 - contig1.name.txt and contig1.fa
	#
	# Input files
	local assembly_graph=${_polap_var_ga_assembly_graph_gfa}
	local assembly_graph_edges_fasta=${_polap_var_ga_assembly_graph_edges_fasta}
	local contig_file=${_polap_var_mtcontigname}
	local output_fasta=${_polap_var_oga_contig}/contig1.fa
	local output_name=${_polap_var_oga_contig}/contig1.name.txt

	>"$output_fasta"

	_polap_log3_cmd mkdir -p ${_polap_var_oga_contig}

	# Step 1: Convert GFA file to FASTA format using gfatools
	gfatools gfa2fa "$assembly_graph" >${assembly_graph_edges_fasta} 2>${_polap_output_dest}

	# Step 2: Loop through each line in mt.contig.name-2 and concatenate the sequences
	line_number=1
	while IFS=, read -r line; do
		# Trim leading/trailing spaces
		line=$(echo $line | xargs)
		_polap_log3 "($line_number): $line"

		# Initialize an empty string to hold the concatenated sequence
		concatenated_sequence=""

		# Process each edge ID in the line
		for edge in $(echo $line | tr ',' '\n'); do
			edge=$(echo $edge | xargs) # Remove extra spaces
			if [[ "$edge" =~ \+$ ]]; then
				_polap_log3 "  (+) $edge"
				# Extract the forward strand sequence (edge_number+)
				edge_id=${edge%+} # Remove the '+' to match the reverse strand
				sequence=$(seqkit grep -p "$edge_id" ${assembly_graph_edges_fasta} | seqkit seq -s)
			elif [[ "$edge" =~ \-$ ]]; then
				_polap_log3 "  (-) $edge"
				# Extract the reverse complementary strand sequence (edge_number-)
				edge_id=${edge%-} # Remove the '-' to match the reverse strand
				sequence=$(seqkit grep -p "$edge_id" ${assembly_graph_edges_fasta} | seqkit seq -v -t DNA -rp | seqkit seq -s)
			fi
			concatenated_sequence+="$sequence"
		done

		# Write the concatenated sequence to the output FASTA file with the line number as the ID
		echo ">contig_${line_number}" >>"$output_fasta"
		echo "$concatenated_sequence" >>"$output_fasta"

		# Increment the line number for the next iteration
		((line_number++))

	done <"$contig_file"

	seqkit seq -ni $output_fasta >$output_name

	# Main 2 - contig2.name.txt and contig2.fa

	# Input files
	local assembly_graph=${_polap_var_ga_assembly_graph_gfa}
	local assembly_graph_edges_fasta=${_polap_var_ga_assembly_graph_edges_fasta}
	local contig_file=${_polap_var_mtcontigname}
	local output_fasta=${_polap_var_oga_contig}/contig2.fa
	local output_name=${_polap_var_oga_contig}/contig2.name.txt
	local temp_fasta=${_polap_var_ga_assembly_graph_edges_fasta}

	# Step 1: Convert GFA file to FASTA format using gfatools
	gfatools gfa2fa "$assembly_graph" >"$temp_fasta" 2>${_polap_output_dest}

	# Step 2: Identify the edge IDs in mt.contig.name-2
	declare -A edge_count
	# Read each line from mt.contig.name-2
	while IFS=, read -r line; do
		# Process each edge ID in the line
		for edge in $(echo $line | tr ',' '\n'); do
			edge=$(echo $edge | xargs) # Remove extra spaces
			if [[ "$edge" =~ \+$ ]]; then
				_polap_log3 "  (+) $edge"
				# Extract the forward strand sequence (edge_number+)
				edge_base=${edge%+} # Remove the '+' to match the reverse strand
			elif [[ "$edge" =~ \-$ ]]; then
				_polap_log3 "  (-) $edge"
				# Extract the reverse complementary strand sequence (edge_number-)
				edge_base=${edge%-} # Remove the '-' to match the reverse strand
			fi
			edge_sign=${edge: -1} # Get the sign (+ or -)

			# Track counts of each edge number
			if [[ -z "${edge_count[$edge_base]+x}" ]]; then
				edge_count[$edge_base]=""
			fi
			edge_count[$edge_base]+="$edge_sign"
		done
	done <"$contig_file"

	# Step 3: Create contig2.name.txt and contig2.fa files
	>"$output_fasta" # Clear the output FASTA file
	>"$output_name"  # Clear the output name file

	# Loop through the edge_count associative array
	for edge_base in "${!edge_count[@]}"; do
		edge_signs="${edge_count[$edge_base]}"

		# Check if both '+' and '-' are present; skip if both are present
		if [[ ! "$edge_signs" =~ "+" ]] || [[ ! "$edge_signs" =~ "-" ]]; then
			# If only one sign is present, select that edge
			if [[ "$edge_signs" =~ "+" ]]; then
				# Select the forward strand
				echo ">${edge_base}+" >>"$output_fasta"
				seqkit grep -p "$edge_base" "$temp_fasta" | seqkit seq -s >>"$output_fasta"
				echo "${edge_base}+" >>"$output_name"
			elif [[ "$edge_signs" =~ "-" ]]; then
				# Select the reverse strand and reverse complement the sequence
				echo ">${edge_base}-" >>"$output_fasta"
				seqkit grep -p "$edge_base" "$temp_fasta" | seqkit seq -v -t DNA -rp -s >>"$output_fasta"
				echo "${edge_base}-" >>"$output_name"
			fi
		fi
	done

	seqkit seq -ni $output_fasta >$output_name

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_map-directional-reads { # selects reads mapped on a genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
This tool maps long reads on seed contigs.

Inputs
------

- source assembly number: i
- destination assembly number: j
- j/01-contig/contig1.fa
- j/01-contig/contig2.fa
- ${_polap_var_outdir_lk_fq_gz}

Outputs
-------

- ${_polap_var_oga_contig}/contig.fa
- ${_polap_var_oga_contig}/contig1.paf
- ${_polap_var_oga_contig}/contig1.tab

Usage
-----
$(basename "$0") ${_arg_menu[0]} -i 1 -j 2
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	if [[ "${_arg_menu[1]}" == "preview" ]]; then
		seqkit stats "${_polap_var_oga_contig}/contig1.fa" >&3
		seqkit stats "${_polap_var_oga_contig}/contig2.fa" >&3

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "view" ]]; then
		ls -l "${_polap_var_oga_contig}/contig1.tab" >&3
		ls -l "${_polap_var_oga_contig}/contig2.tab" >&3

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	_polap_log0 "mapping long-read data on the seed contigs ..."
	_polap_log1 "  assembly: ${_arg_inum} (source) -> ${_arg_jnum} (target) ..."
	_polap_log1 "  input1: ${_polap_var_mtcontigname}"
	_polap_log1 "  input2: ${_polap_var_ga_contigger_edges_fasta}"

	if [[ -s "${_polap_var_oga_contig}/contig1.tab" ]] && [[ "${_arg_redo}" == "on" ]]; then
		_polap_log0 "  found: ${_polap_var_oga_reads}/contig1.tab, so skipping mapping long-read data ..."
		return
	fi

	if [ ! -s "${_polap_var_oga_contig}/contig1.fa" ]; then
		_polap_log0 "ERROR: no such ${_polap_var_oga_contig}/contig1.fa"
		exit $EXIT_ERROR
	fi

	_polap_log1 "  creates ${_polap_var_oga_contig}"
	_polap_log3_cmd mkdir -p "${_polap_var_oga_contig}"

	_polap_log1 "  determines which long-read data to use ..."

	local _source_long_reads_fq=""
	_polap_oga_determine-long-read-file _source_long_reads_fq

	# TODO: 30-contigger/graph_final.gfa must be exist.
	if ! _polap_gfatools-gfa2fasta; then
		_polap_error_message $?
		return ${_POLAP_ERR_MENU_MAP_READS}
	fi

	_polap_log1 "  finding the length of all seed contigs"
	_polap_utility_get_contig_length \
		"${_polap_var_oga_contig}/contig1.fa" \
		"${_polap_var_oga_contig}/contig1_total_length.txt"
	local CONTIG_LENGTH=$(<"${_polap_var_oga_contig}/contig1_total_length.txt")
	_polap_log2_cat "${_polap_var_oga_contig}/contig1_total_length.txt"
	local _contig_length_bp=$(_polap_utility_convert_bp ${CONTIG_LENGTH})
	_polap_log1 "    organelle genome size based on the seed contig selection: ${_contig_length_bp}"

	_polap_log1 "  mapping long-read data on the seed contigs using minimap2 ..."
	_polap_log2 "    input1: ${_polap_var_oga_contig}/contig1.fa"
	_polap_log2 "    input2: ${_source_long_reads_fq}"
	_polap_log2 "    output: ${_polap_var_oga_contig}/contig1.paf"
	if [[ -s "${_polap_var_oga_contig}"/contig1.paf ]] && [[ "${_arg_redo}" = "off" ]]; then
		_polap_log1 "  found: ${_polap_var_oga_reads}/contig1.paf, skipping the minimap2 mapping step ..."
	else
		_polap_log3_pipe "minimap2 -cx \
      ${_arg_minimap2_data_type} \
      ${_polap_var_oga_contig}/contig1.fa \
      ${_source_long_reads_fq} \
      -t ${_arg_threads} \
      -o ${_polap_var_oga_contig}/contig1.paf \
      >${_polap_output_dest} 2>&1"
	fi

	_polap_log1 "  converting PAF to TAB ..."
	_polap_log2 "    input1: ${_polap_var_oga_contig}/contig1.paf"
	_polap_log2 "    output: ${_polap_var_oga_contig}/contig1.tab"
	_polap_log3_cmd bash "$script_dir/run-polap-sh-minimap2-paf2tab.sh" "${_arg_min_read_length}" "${_polap_var_oga_contig}/contig1.paf" "${_polap_var_oga_contig}/contig1.tab"

	_polap_log1 "  mapping long-read data on the seed contigs using minimap2 ..."
	_polap_log2 "    input1: ${_polap_var_oga_contig}/contig2.fa"
	_polap_log2 "    input2: ${_source_long_reads_fq}"
	_polap_log2 "    output: ${_polap_var_oga_contig}/contig2.paf"
	if [[ -s "${_polap_var_oga_contig}"/contig2.paf ]] && [[ "${_arg_redo}" = "off" ]]; then
		_polap_log1 "  found: ${_polap_var_oga_reads}/contig2.paf, skipping the minimap2 mapping step ..."
	else
		_polap_log3_pipe "minimap2 -cx \
      ${_arg_minimap2_data_type} \
      ${_polap_var_oga_contig}/contig2.fa \
      ${_source_long_reads_fq} \
      -t ${_arg_threads} \
      -o ${_polap_var_oga_contig}/contig2.paf \
      >${_polap_output_dest} 2>&1"
	fi

	_polap_log1 "  converting PAF to TAB ..."
	_polap_log2 "    input1: ${_polap_var_oga_contig}/contig2.paf"
	_polap_log2 "    output: ${_polap_var_oga_contig}/contig2.tab"
	_polap_log3_cmd bash "$script_dir/run-polap-sh-minimap2-paf2tab.sh" "${_arg_min_read_length}" "${_polap_var_oga_contig}/contig2.paf" "${_polap_var_oga_contig}/contig2.tab"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
#
################################################################################
function _run_polap_test-directional-reads { # selects reads mapped on a genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
This tool tests organelle-genome assemblies with directional reads.

2. select read names: 02-reads
3. select reads: 03-seeds
4. subsample reads: 04-subsample
5. flye run: 05-flye
6. summary: 06-summary
7. plot: 07-plot

Inputs
------

- source assembly number: i
- destination assembly number: j
- -s, --select-read-range <start,end,count>
- --start-index <index>: to start at somewhere not start
- input long read data to use (priority order): 
  1. ${_polap_var_outdir_lk_fq_gz}
  2. ${_polap_var_outdir_nk_fq_gz}
  3. ${_arg_long_reads}

Outputs
-------

- ${_polap_var_oga_contig}: map-reads output
- ${_polap_var_oga_reads}: mapped read names
- ${_polap_var_oga_seeds}: mapped reads
- ${_polap_var_oga_subsample}: subsample of the mapped reads
- ${_polap_var_oga_flye}: assemblies
- ${_polap_var_oga_summary}: summary
- ${_polap_var_oga_plot}: plot or table using range in contig folder and summary

Submenu
-------

- ptgaul
- polap

View
----

- ptgaul
- polap

Report
------

- ptgaul
- polap

Usage
-----
$(basename "$0") ${_arg_menu[0]} [ptgaul] --select-read-range 3000,27000,5
$(basename "$0") ${_arg_menu[0]} polap -i 1 -j 2
$(basename "$0") ${_arg_menu[0]} polap --select-read-range 3000,27000,5
$(basename "$0") ${_arg_menu[0]} polap --select-read-range 3000,27000,5 --start-index 3
$(basename "$0") ${_arg_menu[0]} view polap -i 2 -j 3
$(basename "$0") ${_arg_menu[0]} report ptgaul --report-x 3000,5000,7000,9000,11000
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	if [[ "${_arg_menu[1]}" == "preview" ]]; then
		ls -l "${_polap_var_oga_contig}/contig1.tab" >&3
		ls -l "${_polap_var_oga_contig}/contig2.tab" >&3

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "view" ]]; then
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "report" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	_polap_log0 "test directional read-selection ..."

	##############################################################################
	# Set the default submenu option for the main menu
	#   _pread_sel: submenu
	#   _read_names: read name file to use (ptgaul, single, combined)
	if [[ "${_arg_menu[1]}" == "infile" ]]; then
		_arg_menu[1]="ptgaul"
		_polap_log0 "  default set to read-selection: ${_arg_menu[1]}"
	fi
	local _pread_sel="${_arg_menu[1]}"
	local _read_names="${_arg_menu[1]}"

	_polap_log1 "  argument1: ${_pread_sel}"
	_polap_log1 "  argument2: ${_read_names}"

	##############################################################################
	# Main input files: seed contigs and long-read data to map on the contigs
	#   The Minimap2 mapping should have been completed
	_polap_log1 "  input1: ${_polap_var_oga_contig}"

	check_folder_existence "${_polap_var_oga_contig}"

	_polap_utility_get_contig_length \
		"${_polap_var_oga_contig}/contig1.fa" \
		"${_polap_var_oga_contig}/contig1_total_length.txt"
	local CONTIG_LENGTH=$(<"${_polap_var_oga_contig}/contig1_total_length.txt")
	_polap_log2_cat "${_polap_var_oga_contig}/contig1_total_length.txt"

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
	_create_range "${_arg_select_read_range}" \
		"${_polap_var_oga_contig}/${_pread_sel}.txt"
	_polap_log2_cat "${_polap_var_oga_contig}/${_pread_sel}.txt"

	# Check if one wants to replace the previous run
	# if [[ -d "${_polap_var_oga_reads}/${_pread_sel}" ]] &&
	# 	[[ "${_arg_start_index}" -eq 0 ]] &&
	# 	[[ "${_arg_redo}" == "off" ]]; then
	# 	if confirm "Do you want to redo the test-reads of ${_pread_sel} on assembly ${_arg_jnum}"?; then
	# 		_arg_redo="on"
	# 	else
	# 		_polap_log0 "You cancelled the operation."
	# 		return 0
	# 	fi
	# fi

	if [[ "${_arg_redo}" == "on" ]]; then
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

		_arg_single_min="${restored_array[i]}"
		_arg_pair_min="${restored_array[i]}"
		_arg_bridge_min="0"

		_polap_log3_cmd mkdir -p "${_polap_var_oga_reads}/${_pread_sel}/${i}"

		if [[ -s "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz" ]]; then
			_polap_log0 "  found: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz, so skipping ..."
		else
			_polap_log0 "  (${i}) selecting read names for ${_pread_sel}: w = ${_test_value} (bp)"
			_polap_log2 "    single-min: ${_arg_single_min}"
			_polap_log2 "    pair-min: ${_arg_pair_min}"
			_polap_log2 "    bridge-min: ${_arg_bridge_min}"
			_polap_log2 "    input1: ${_polap_var_oga_contig}/contig1.tab"
			_polap_log2 "    input2: ${_polap_var_oga_contig}/contig2.tab"
			_polap_log2 "    output1: ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.contig1.forward.names"
			_polap_log2 "    output1: ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.contig2.forward.names"
			_polap_log2 "    output1: ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.contig.forward.names"
			_polap_log2 "    output1: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.forward.fq.gz"
			_polap_log2 "    output1: ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.contig2.reverse.names"
			_polap_log2 "    output1: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.reverse.fq.gz"
			_polap_log2 "    output1: ${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"

			# contig1 -> same strand
			# contig2 -> reverse complementary

			_polap_log3_pipe "Rscript --vanilla ${script_dir}/run-polap-r-select-reads-ptgaul.R \
        --forward \
		    -t ${_polap_var_oga_contig}/contig1.tab \
		    -o ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.contig1.forward.names.txt \
        -w ${_arg_single_min} \
		    >${_polap_output_dest} 2>&1"

			# _polap_log3_pipe "Rscript --vanilla ${script_dir}/run-polap-r-select-reads-ptgaul.R \
			#      --forward \
			#    -t ${_polap_var_oga_contig}/contig2.tab \
			#    -o ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.contig2.forward.names.txt \
			#      -w ${_arg_single_min} \
			#    >${_polap_output_dest} 2>&1"
			#
			_polap_log3_pipe "cat \
			  ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.contig1.forward.names.txt |\
			  sort | uniq \
			  >${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.contig.forward.names.txt"

			_polap_log3_pipe "seqtk subseq \
			   ${_source_long_reads_fq} \
			   ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.contig.forward.names.txt |\
			   gzip >${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"

			# _polap_log3_pipe "Rscript --vanilla ${script_dir}/run-polap-r-select-reads-ptgaul.R \
			#      --reverse \
			#    -t ${_polap_var_oga_contig}/contig2.tab \
			#    -o ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.contig2.reverse.names.txt \
			#      -w ${_arg_single_min} \
			#    >${_polap_output_dest} 2>&1"
			#
			# # delete reverse names that belong to the forward name group.
			# _polap_log3_pipe "grep -v \
			# 	-f ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.contig1.forward.names.txt \
			# 	${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.contig2.reverse.names.txt \
			# 	>${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.contig2.reverse.proper.names.txt"
			#
			# _polap_log3_pipe "seqtk subseq \
			#    ${_source_long_reads_fq} \
			#    ${_polap_var_oga_reads}/${_pread_sel}/${i}/${_read_names}.contig2.reverse.proper.names.txt |\
			#      seqtk seq -r |\
			#    gzip >${_polap_var_oga_seeds}/${_pread_sel}/${i}.reverse.fq.gz"
			#
			# _polap_log3_pipe "zcat \
			#   ${_polap_var_oga_seeds}/${_pread_sel}/${i}.forward.fq.gz \
			#   ${_polap_var_oga_seeds}/${_pread_sel}/${i}.reverse.fq.gz |
			#   gzip >${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"
		fi
	done

	# Iterate over the array using an index
	for ((i = ${_arg_start_index}; i < array_length; i++)); do
		local _test_value="${restored_array[i]}"

		if [[ -s "${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz" ]]; then
			_polap_log0 "  found: ${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz, so skipping ..."
		else
			# subsample
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
            -o ${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz \
            2>${_polap_output_dest}"
					_polap_log3_pipe "echo ${_random_seed} >${_polap_var_oga_subsample}/${_pread_sel}/${i}.random.seed.${_random_seed}"
				else
					_polap_log0 "    no reduction of the long-read data because of the option --no-coverage-check: expected coverage: ${_expected_organelle_coverage}x"
					_polap_log3_cmd ln -s $(realpath "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz") "${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz"
				fi
			else
				_polap_log0 "    no reduction of the long-read data because $_expected_organelle_coverage < ${_arg_coverage}"
				_polap_log3_cmd ln -s $(realpath "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz") "${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz"
			fi
		fi
	done

	return 0

	# Iterate over the array using an index
	for ((i = ${_arg_start_index}; i < array_length; i++)); do
		local _test_value="${restored_array[i]}"

		if [[ -s "${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa" ]]; then
			_polap_log0 "  found: ${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa, so skipping ..."
		else
			_polap_log1 "  dflye assembly for ${_pread_sel}"
			_polap_log2 "    input1: ${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz"
			_polap_log2 "    output: ${_polap_var_oga_flye}/${_pread_sel}/${i}/30-contigger/graph_final.gfa"
			local _command1="dflye \
        ${_arg_flye_data_type} \
        ${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz \
		    --out-dir ${_polap_var_oga_flye}/${_pread_sel}/${i} \
        --directional-reads \
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
		fi
	done

	return 0

	# Read the file contents into an array
	read -a restored_array <"${_polap_var_oga_contig}/${_pread_sel}.txt"
	local array_length=${#restored_array[@]}
	# Iterate over the array using an index
	for ((i = ${_arg_start_index}; i < array_length; i++)); do
		local _test_value="${restored_array[i]}"

		_arg_single_min="${restored_array[i]}"
		_arg_pair_min="${restored_array[i]}"
		_arg_bridge_min="0"

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

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_select-directional-reads { # selects reads mapped on a genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<HEREDOC
Test organelle-assembly results with a chosen w value multiple times.


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

$(basename "$0") ${_arg_menu[0]} [ptgaul-reads] -w 3000
$(basename "$0") ${_arg_menu[0]} intra-reads -w 5000
$(basename "$0") ${_arg_menu[0]} polap-reads -w 3000 5
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
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_original-directional-reads { # selects reads mapped on a genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	source "$script_dir/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
This tool tests directional reads.

2. select read names: 02-reads
3. select reads: 03-seeds
4. subsample reads: 04-subsample
5. flye run: 05-flye
6. summary: 06-summary
7. plot: 07-plot

# Map long reads on a Flye genome assembly.
#
# Arguments:
#   -l ${_arg_long_reads}: long-read data default:${_polap_var_outdir_lk_fq_gz}
#   -p ${_arg_unpolished_fasta}
#   -f ${_arg_final_assembly}
# Inputs:
#   ${_polap_var_ga_contigger_edges_fasta}
#   ${_polap_var_outdir_lk_fq_gz} or ${_polap_var_outdir_nk_fq_gz}
# Outputs:
#   ${_polap_var_oga_contig}/seed.fq
Example: $(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads} -p seed.fa -f seed.fq

HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	mkdir -p ${_polap_var_oga_contig}
	_polap_log3_cmd cp ${_arg_unpolished_fasta} ${_polap_var_oga_contig}/contig.fa
	local _source_long_reads_fq=""
	_polap_oga_determine-long-read-file _source_long_reads_fq
	_polap_log1 "  finding the length of all seed contigs"
	_polap_utility_get_contig_length \
		"${_polap_var_oga_contig}/contig.fa" \
		"${_polap_var_oga_contig}/contig_total_length.txt"
	local CONTIG_LENGTH=$(<"${_polap_var_oga_contig}/contig_total_length.txt")

	_pread_sel="directional"
	i="1"

	if [[ "${_arg_menu[1]}" == "1" ]]; then

		_polap_log0 "mapping long-read data on the seed contigs ..."

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
		_polap_log3_cmd bash "$script_dir/run-polap-sh-minimap2-paf2tab.sh" "${_arg_min_read_length}" "${_polap_var_oga_contig}/contig.paf" "${_polap_var_oga_contig}/contig.tab"

	fi

	if [[ "${_arg_menu[1]}" == "2" ]]; then
		mkdir -p ${_polap_var_oga_reads}/${_pread_sel}/${i}
		mkdir -p ${_polap_var_oga_seeds}/${_pread_sel}/${i}
		mkdir -p ${_polap_var_oga_subsample}/${_pread_sel}/${i}

		_polap_log3_pipe "Rscript --vanilla ${script_dir}/run-polap-r-directional.R \
        --use-strand \
	    -m ${_polap_var_mtcontigname} \
		  -t ${_polap_var_oga_contig}/contig.tab \
		  --out ${_polap_var_oga_reads}/${_pread_sel}/${i} \
      -w ${_arg_single_min} \
		  -r ${_arg_pair_min} \
		  -x ${_arg_bridge_min} \
      --create-ptgaul \
		  >${_polap_output_dest} 2>&1"

		_polap_log3_pipe "seqtk subseq \
		    ${_source_long_reads_fq} \
		    ${_polap_var_oga_reads}/${_pread_sel}/${i}/ptgaul.names |\
		    gzip >${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz"
	fi

	if [[ "${_arg_menu[1]}" == "3" ]]; then

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
        -o ${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz \
        2>${_polap_output_dest}"
				_polap_log3_pipe "echo ${_random_seed} >${_polap_var_oga_subsample}/${_pread_sel}/${i}.random.seed.${_random_seed}"
			else
				_polap_log0 "    no reduction of the long-read data because of the option --no-coverage-check: expected coverage: ${_expected_organelle_coverage}"
				_polap_log3_cmd ln -s $(realpath "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz") "${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz"
			fi
		else
			_polap_log0 "    no reduction of the long-read data because $_expected_organelle_coverage < ${_arg_coverage}"
			_polap_log3_cmd ln -s $(realpath "${_polap_var_oga_seeds}/${_pread_sel}/${i}.fq.gz") "${_polap_var_oga_subsample}/${_pread_sel}/${i}.fq.gz"
		fi
	fi

	_polap_log1 "NEXT: $(basename "$0") reads -o ${_arg_outdir} -i ${_arg_inum} -j ${_arg_jnum}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
