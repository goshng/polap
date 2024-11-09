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

# depth_range=()
# _polap_seeds_get-depth-range-of depth-range.txt depth_range
# ${depth_raneg[0]}
# ${depth_raneg[1]}
function _polap_seeds_get-depth-range-of() {
	local _depth_range=$1
	local -n _r_arr=$2
	# Check if there is a manual depth range file.
	if [[ -s "${_depth_range}" ]]; then
		# Extract numbers from the file
		local numbers=$(awk 'NR==2 {print $1, $2}' "${_depth_range}")
		# Store the numbers in variables
		local depth_lower=$(echo $numbers | cut -d ' ' -f 1)
		local depth_upper=$(echo $numbers | cut -d ' ' -f 2)
		_polap_log2 "    custom depth range: $depth_lower ~ $depth_upper"
		_polap_log3_column "${_depth_range}"
	else
		_polap_log0 "ERROR: no such file: ${_depth_range}"
		local depth_lower=0
		local depth_upper=0
	fi
	_r_arr=($depth_lower $depth_upper)
}

################################################################################
# Create a depth range file based on user input.
# If a user supplies a text file containing two numerical values, then generate
# a depth range file accordingly.
# If not using an existing depth range file, then request that users provide
# two specific values to generate such a file.
################################################################################
function _polap_seeds_create-manual-depth-range() {
	local _two_value_textfile="${_arg_menu[2]}"

	if [[ -s "${_two_value_textfile}" ]]; then
		_polap_log1 "    input3 (two-value file): ${_two_value_textfile}"
		_polap_log2_cat "${_two_value_textfile}"
		readarray -t _numbers < <(grep -Eo '[0-9]+' "${_two_value_textfile}")
		_polap_log1 "    manual depth range: ${_numbers[0]}"
		_polap_log1 "    manual depth range: ${_numbers[0]} ~ ${_numbers[1]}"
		_polap_log3_cmd bash "$script_dir/run-polap-sh-create-depth-file.sh" "${_numbers[0]}" "${_numbers[1]}" "${_polap_var_mtcontigs_1_custom_depth_range}"
		_polap_log2_column "${_polap_var_mtcontigs_1_custom_depth_range}"
		_polap_log1 "    use this depth-range for both contig preselection and graph filtering"
		_polap_log3_cmd cp "${_polap_var_mtcontigs_1_custom_depth_range}" "${_polap_var_mtcontigs_2_custom_depth_range}"
	else
		# 1. for the depth-range in the contig preselection
		_polap_log0 "  Choose desired values of the depth range for the contig preselection:"
		_polap_log0 "    we will use the depth-range for the graph filtering if you do not have another depth-range file."
		_polap_log0 "    1. First, calculate an estimated average value accordingly."
		_polap_log0 "    2. Choosing a lower bound that is one-third of the average value can be a suitable option."
		_polap_log0 "    3. Selecting an appropriate upper limit can be achieved by setting it to three times the average."
		read -p "  Enter the lower bound of a depth range: " _num1
		read -p "  Enter the upper bound of a depth range: " _num2
		_polap_log3_cmd bash "$script_dir/run-polap-sh-create-depth-file.sh" "${_num1}" "${_num2}" "${_polap_var_mtcontigs_1_custom_depth_range}"
		_polap_log0 "    manual depth range: ${_num1} ~ ${_num2}"
		_polap_log2_column "${_polap_var_mtcontigs_1_custom_depth_range}"
		# 2. for the depth-range in the graph filtering
		if [[ "${_two_value_textfile}" == "2" ]]; then
			# Prompt the user to enter the first number
			_polap_log0 "  Choose desired values of the depth range for the graph filtering:"
			_polap_log0 "    1. First, calculate an estimated average value accordingly."
			_polap_log0 "    2. Choosing a lower bound that is one-third of the average value can be a suitable option."
			_polap_log0 "    3. Selecting an appropriate upper limit can be achieved by setting it to three times the average."
			read -p "  Enter the lower bound of a depth range: " _num1
			read -p "  Enter the upper bound of a depth range: " _num2
			_polap_log3_cmd bash "$script_dir/run-polap-sh-create-depth-file.sh" "${_num1}" "${_num2}" "${_polap_var_mtcontigs_2_custom_depth_range}"
			_polap_log0 "    manual depth range: ${_num1} ~ ${_num2}"
			_polap_log2_column "${_polap_var_mtcontigs_2_custom_depth_range}"
		else
			_polap_log1 "    use this depth-range for both contig preselection and graph filtering"
			_polap_log3_cmd cp "${_polap_var_mtcontigs_1_custom_depth_range}" "${_polap_var_mtcontigs_2_custom_depth_range}"
		fi
	fi
}

################################################################################
#
# determines the depth range automatically ..."
# using the length distribution of contigs sorted"
# by copy numbers and organelle gene counts"
################################################################################
function _polap_seeds_create-automatic-depth-range() {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	local _type=$1
	local -n result_ref=$2

	if [[ "${_arg_plastid}" == "on" ]]; then
		_polap_log2 "  using plastid depth range ..."
		local _command1="Rscript $script_dir/run-polap-r-plastid-determine-depth-range_${_type}.R \
				-t ${_polap_var_annotation_table} \
				-c ${_polap_var_ga_annotation_cdf_table} \
				-o ${_polap_var_mtcontigs_1_custom_depth_range} \
        --plastid \
				2>$_polap_output_dest"
	else
		local _command1="Rscript $script_dir/run-polap-r-determine-depth-range_${_type}.R \
				-t ${_polap_var_annotation_table} \
				-c ${_polap_var_ga_annotation_cdf_table} \
				-o ${_polap_var_mtcontigs_1_custom_depth_range} \
				2>$_polap_output_dest"
	fi
	_polap_log3_pipe "${_command1}"

	_polap_log3_column "${_polap_var_ga_annotation_cdf_table}"

	if [[ -s "${_polap_var_mtcontigs_1_custom_depth_range}" ]]; then
		_polap_log1 "    we use one depth-range for contig preselection and graph filtering"
		_polap_log3_pipe "cp ${_polap_var_mtcontigs_1_custom_depth_range} \
      ${_polap_var_mtcontigs_2_custom_depth_range}"
		_polap_log2 "    Automatically generated depth range!"
		result_ref="depth-range"
	else
		_polap_log2 "    No automatically generated depth range!"
		result_ref="no depth-range"
	fi
}

################################################################################
# Pre-select contigs that have been identified to originate from the
# mitochondrial genome.
################################################################################
function _polap_seeds_preselect-contigs() {

	if [[ -s "${_polap_var_mtcontigs_depth_range_preselection}" ]]; then
		# Extract numbers from the file
		local numbers=$(awk 'NR==2 {print $1, $2}' "${_polap_var_mtcontigs_depth_range_preselection}")
		# Store the numbers in variables
		local depth_lower=$(echo $numbers | cut -d ' ' -f 1)
		local depth_upper=$(echo $numbers | cut -d ' ' -f 2)
		_polap_log0 "  depth range for edge contig preselection from ${_polap_var_mtcontigs_depth_range_preselection}: $depth_lower ~ $depth_upper"
		_polap_log2_column "${_polap_var_mtcontigs_depth_range_preselection}"
	else
		die "ERROR: no such file: ${_polap_var_mtcontigs_depth_range_preselection}"
	fi

	_polap_log2 "    minimum gene density for mtDNA: 10 per 1 Mb"
	_polap_log2 "    gene count comparison: MT > PT"
	_polap_log2 "    depth range: $depth_lower ~ $depth_upper"
	_polap_log2 "    input1: ${_polap_var_mtcontigs_depth_range_preselection}"
	_polap_log2 "    input2: ${_polap_var_annotation_table}"
	_polap_log2 "    output: ${_polap_var_mtcontigs_preselection}"
	local _command1="Rscript $script_dir/run-polap-r-preselect-annotation.R  \
		--out ${_polap_var_mtcontigs_preselection} \
		--table ${_polap_var_annotation_table} \
		--depth-range ${depth_lower},${depth_upper} \
		--compare-mt-pt \
		--gene-density 10"
	if [[ "${_arg_plastid}" = "on" ]]; then
		_polap_log2 "    for plastid not mitochondrial DNA"
		_command1+=" \
        --plastid"
	fi
	_command1+=" \
				2>$_polap_output_dest"
	_polap_log3_pipe "${_command1}"

	if [[ -s "${_polap_var_mtcontigs_preselection}" ]]; then
		_polap_log2_column "${_polap_var_mtcontigs_preselection}"
	else
		_polap_log0 "ERROR: no such file: ${_polap_var_mtcontigs_preselection}"
	fi
}

################################################################################
################################################################################
function _polap_seeds_depthfilter-gfa() {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	local depth_values=()
	_polap_seeds_get-depth-range-of \
		"${_polap_var_mtcontigs_depth_range_graphfilter}" \
		depth_values

	if [[ "${depth_values[0]}" -gt 0 ]]; then
		local depth_lower="${depth_values[0]}"
		local depth_upper="${depth_values[1]}"
		_polap_log2 "  depth range for graph filtering: $depth_lower ~ $depth_upper"
	else
		die "ERROR: no depth ranges"
	fi

	_polap_log1 "    step 3-1: creating GFA without sequence data: ${_polap_var_mtcontigs_gfa_all}"
	_polap_log1 "      input1: ${_polap_var_ga_contigger_edges_gfa}"
	_polap_log2 "      output: ${_polap_var_mtcontigs_gfa_all}"
	_polap_log3_pipe "gfatools view \
		-S ${_polap_var_assembly_graph_final_gfa} \
		>${_polap_var_mtcontigs_gfa_all} \
		2>$_polap_output_dest"

	_polap_log1 "    step 3-2: extracting sequence part of GFA: ${_polap_var_mtcontigs_gfa_seq_part}"
	_polap_log2 "      input1: ${_polap_var_mtcontigs_gfa_all}"
	_polap_log2 "      output: ${_polap_var_mtcontigs_gfa_seq_part}"
	_polap_log3_pipe "grep ^S ${_polap_var_mtcontigs_gfa_all} >${_polap_var_mtcontigs_gfa_seq_part}"

	# Filter edges in GFA using depths.
	_polap_log1 "    step 3-3: filtering GFA sequence part using depth range"
	_polap_log2 "      input1: ${_polap_var_mtcontigs_gfa_seq_part}"
	_polap_log2 "      input2: ${_polap_var_mtcontigs_depth_range_graphfilter}"
	_polap_log2 "        depth range: $depth_lower ~ $depth_upper"
	_polap_log2 "      output: ${_polap_var_mtcontigs_gfa_seq_filtered}"
	if [[ -s "${_polap_var_mtcontigs_depth_range_graphfilter}" ]]; then
		_polap_log3_pipe "Rscript $script_dir/run-polap-r-depthfilter-gfa.R \
			--gfa ${_polap_var_mtcontigs_gfa_seq_part} \
			--depth ${_polap_var_mtcontigs_depth_range_graphfilter} \
			--out ${_polap_var_mtcontigs_gfa_seq_filtered} \
			2>$_polap_output_dest"
	else
		_polap_log0 "    no such file: ${_polap_var_mtcontigs_depth_range_graphfilter}"
		_polap_log0 "    ${MTCONTIGNAME} -> empty or no seed contigs"
		>"${MTCONTIGNAME}"
		_polap_log2_cat "${MTCONTIGNAME}"
		_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	if [[ -s "${_polap_var_mtcontigs_gfa_seq_filtered}" ]]; then
		_polap_log3_head "${_polap_var_mtcontigs_gfa_seq_filtered}"
	fi
}

################################################################################
# Prepare data
################################################################################
function _polap_seeds_prepare-cc() {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	_polap_log2 "    step 4-1: subsetting GFA using the depth-filtered GFA sequence part: ${_polap_var_mtcontigs_gfa_depthfiltered_gfa}"
	_polap_log2 "      input1: ${_polap_var_mtcontigs_gfa_seq_filtered}"
	_polap_log2 "      output: ${_polap_var_mtcontigs_gfa_depthfiltered_gfa}"

	_polap_log3_pipe "cut -f1 \
    ${_polap_var_mtcontigs_gfa_seq_filtered} \
    >${_polap_var_mtcontigs_gfa_seq_filtered_edge}"

	_polap_log3_pipe "gfatools view -S \
		-l @${_polap_var_mtcontigs_gfa_seq_filtered_edge} \
		${_polap_var_assembly_graph_final_gfa} \
		2>$_polap_output_dest \
		>${_polap_var_mtcontigs_gfa_depthfiltered_gfa}"

	if [[ -s "${_polap_var_mtcontigs_gfa_depthfiltered_gfa}" ]]; then
		_polap_log2_head "${_polap_var_mtcontigs_gfa_depthfiltered_gfa}"
	fi

	# Prepare links for finding connected components.
	_polap_log2 "    step 4-2: preparing links for finding connected components"
	_polap_log2 "      input: ${_polap_var_mtcontigs_gfa_depthfiltered_gfa}"
	_polap_log2 "      output: ${_polap_var_mtcontigs_links_tsv}"
	_polap_log3_pipe "grep ^L ${_polap_var_mtcontigs_gfa_depthfiltered_gfa} |\
    cut -f2,4 >${_polap_var_mtcontigs_links_tsv}"

	# No links case
	if [[ -s "${_polap_var_mtcontigs_links_tsv}" ]]; then
		_polap_log2_head "${_polap_var_mtcontigs_links_tsv}"
	else
		# no links -> just current selection
		_polap_log3_pipe "cut -f2 ${_polap_var_mtcontigs_gfa_depthfiltered_gfa} |
			sort | uniq >${MTCONTIGNAME}"
		_polap_log0 "  output1: ${MTCONTIGNAME}"
		_polap_log2_cat "${MTCONTIGNAME}"

		_polap_seeds_final-mtcontig

		_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	local _preselection="${_polap_var_mtcontigs_preselection}"

	# Run R script to analyze GFA links
	_polap_log2 "    step 4-3: preparing for finding connected components"
	_polap_log2 "      input1: ${_preselection}"
	_polap_log2 "      input2: ${_polap_var_mtcontigs_links_tsv}"
	_polap_log2 "      output-base: ${_polap_var_mtcontigs_links}"
	_polap_log2 "      output1: ${_polap_var_mtcontigs_links_number}"
	_polap_log2 "      output2: ${_polap_var_mtcontigs_links_order}"
	_polap_log2 "      output3: ${_polap_var_mtcontigs_links_contig}"
	_polap_log2 "      output4: ${_polap_var_mtcontigs_links_contig_na}"
	_polap_log3_pipe "Rscript $script_dir/run-polap-r-prepare-cc.R \
		--edge ${_preselection} \
		--gfa ${_polap_var_mtcontigs_links_tsv} \
		--out ${_polap_var_mtcontigs_links} \
		--unit-annotation-edge \
		2>$_polap_output_dest"

	if [[ -s "${_polap_var_mtcontigs_links_number}" ]]; then
		_polap_log2_head "${_polap_var_mtcontigs_links_number}"
	fi
	if [[ -s "${_polap_var_mtcontigs_links_order}" ]]; then
		_polap_log2_head "${_polap_var_mtcontigs_links_order}"
	fi
	if [[ -s "${_polap_var_mtcontigs_links_contig}" ]]; then
		_polap_log2_head "${_polap_var_mtcontigs_links_contig}"
	fi
	if [[ -s "${_polap_var_mtcontigs_links_contig_na}" ]]; then
		_polap_log2_head "${_polap_var_mtcontigs_links_contig_na}"
	else
		_polap_log2 "No NA edge: ${_polap_var_mtcontigs_links_contig_na}"
	fi
}

################################################################################
# Finalize the mt.contig.name.
################################################################################
function _polap_seeds_final-mtcontig() {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# local _is_auto=$1 not using it anymore? for --length 1000000 option below

	_polap_log1 "  step 8-1: filtering by max length 1 Mb of contig seeds"
	_polap_log1 "    1 Mb length filtering for auto"
	_polap_log2 "    input1: ${_polap_var_annotation_table}"
	_polap_log2 "    input2: ${_polap_var_mtcontigs_7mtcontigname}"
	_polap_log2 "    output1: ${_polap_var_mtcontig_table}"
	local _command1="Rscript $script_dir/run-polap-r-final-filter-mtcontig.R \
		-t ${_polap_var_annotation_table} \
		-m ${_polap_var_mtcontigs_7mtcontigname} \
    -o ${_polap_var_mtcontig_table} \
    --length 1000000 \
		2>$_polap_output_dest"
	_polap_log3_pipe "${_command1}"

	_polap_log1 "  step 8-2: mt.contig.table with contig seed marks"
	_polap_log2 "    input1: ${_polap_var_annotation_table}"
	_polap_log2 "    input2: ${_polap_var_ga_annotation_depth_table}"
	_polap_log2 "    input3: ${_polap_var_mtcontigs_7mtcontigname}"
	_polap_log2 "    output1: ${_polap_var_mtcontigs_annotation_table_seed}"
	local _command1="Rscript $script_dir/run-polap-r-final-seed-mtcontig.R \
		-t ${_polap_var_annotation_table} \
		-a ${_polap_var_ga_annotation_depth_table} \
		-m ${_polap_var_mtcontigs_7mtcontigname} \
    -o ${_polap_var_mtcontigs_annotation_table_seed} \
		2>$_polap_output_dest"
	_polap_log3_pipe "${_command1}"

}

function _polap_seeds_report-mtcontig() {
	if [[ -s "${_polap_var_mtcontigs_8mtcontigname}" ]]; then
		_polap_log1_cat "${_polap_var_mtcontigs_8mtcontigname}"
		_polap_log0 "---"
		if [[ "${_arg_log_stderr}" = "off" ]]; then
			paste -sd',' "${_polap_var_mtcontigs_8mtcontigname}" >&3
		else
			paste -sd',' "${_polap_var_mtcontigs_8mtcontigname}" >&2
		fi
	fi
}

function _run_polap_choose-seed() { # select seed contigs
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh"
	source "$script_dir/polap-variables-mtcontigs.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Select contigs based on gene density, a customized depth of coverage, and 
# an analysis of the assembly graph.
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   -k $KNUM: destination Flye organelle assembly number
#   --plastid
# Inputs:
#   ${_polap_var_mtcontigs_7mtcontigname}
# Outputs:
#   ${_polap_var_mtcontigname}
Example: $0 ${_arg_menu[0]} -i 1 -j 2 -k 1
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		return
	fi

	# We initiate the process of selecting seed contigs.
	_polap_log0 "choose seed contigs using the assembly graph: ${INUM} (source) -> ${JNUM} (target) ..."

	if [[ -s "${_polap_var_mtcontig_table}" ]]; then
		_polap_log3_column "${_polap_var_mtcontig_table}"
		_polap_log2_pipe "cut -f1 ${_polap_var_mtcontig_table} \
      >${_polap_var_mtcontigname}"
	else
		die "ERROR: no final mtcontig: ${_polap_var_mtcontig_table}"
	fi

	if [[ -s "${_polap_var_mtcontigname}" ]]; then
		_polap_log1_cat "${_polap_var_mtcontigname}"
		_polap_log0 "---"
		if [[ "${_arg_log_stderr}" = "off" ]]; then
			paste -sd',' "${_polap_var_mtcontigname}" >&3
		else
			paste -sd',' "${_polap_var_mtcontigname}" >&2
		fi
	fi

	_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
#
################################################################################
function _run_polap_seeds() { # select seed contigs
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh"
	source "$script_dir/polap-variables-mtcontigs.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Select contigs using seeds-graph, seeds-gene, ...
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   --plastid
# Inputs:
#   ${_polap_var_assembly_graph_final_gfa}
#   ${_polap_var_annotation_table}
# Outputs:
#   ${_polap_var_mtcontigs_8mtcontigname}
# View:
#   <number> for the mt.contig.name-<number>
Example: $0 ${_arg_menu[0]} -i 0 -j 1
Example: $0 ${_arg_menu[0]} view 1 
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		local index="${_arg_menu[2]}"
		local new_filename_table="${_polap_var_ga}/mtcontig-annotation-table-seed-${index}.txt"
		_polap_log0_column "${new_filename_table}"

		local _mt_contig_name="${_polap_var_ga}/mt.contig.name-${index}"
		if [[ -s "${_mt_contig_name}" ]]; then
			_polap_log0 "---"
			if [[ "${_arg_log_stderr}" = "off" ]]; then
				paste -sd',' "${_mt_contig_name}" >&3
			else
				paste -sd',' "${_mt_contig_name}" >&2
			fi
		fi
		return 0
	fi

	# Backup mt.contig.name files.
	if [[ "${_arg_menu[1]}" == "backup" ]]; then
		# Format the current date as YYYY-MM-DD
		local folder_name=$(date +"%Y-%m-%d")
		# Create the folder with the formatted date as its name
		mkdir "${_polap_var_ga}/backup-mt.contig.name-$folder_name"
		mv "${_polap_var_ga}"/mt.contig.name-* "${_polap_var_ga}/backup-mt.contig.name-$folder_name"
		return 0
	fi

	# We initiate the process of selecting seed contigs.
	_polap_log0 "selecting seed contigs using the assembly graph: ${INUM} (source) -> ${JNUM} (target) ..."

	if [[ "${_arg_menu[1]}" = "bandage" ]]; then
		# Prompt the user for input
		read -p "Enter edges using Bandage (e.g., edge_265, edge_520, edge_425): " edges

		# Convert the input string into new lines (replace ", " with newlines)
		local formatted_edges=$(echo "$edges" | tr ', ' '\n' | sed '/^ *$/d')
		echo "${formatted_edges}" >"${MTCONTIGNAME}"
		# _polap_seeds_final-mtcontig
		return $RETURN_SUCCESS
	fi

	_polap_log1 "  input1: ${_polap_var_ga_contigger_edges_gfa}"
	_polap_log1 "  input2: ${_polap_var_annotation_table}"

	# Create the array based on the value of _arg_plastid
	if [[ "$_arg_plastid" == "on" ]]; then
		knum_array=(1 6)
	else
		knum_array=({1..6}) # This creates an array of 1 through 6
	fi

	# Loop over the array and assign each element to KNUM, then call _run_polap_seeds
	for i in "${knum_array[@]}"; do
		KNUM=$i
		_polap_log2 "  running _run_polap_seeds with KNUM=$KNUM"

		# Grouped file path declarations
		source "$script_dir/polap-variables-mtcontigs.sh"
		_run_polap_seeds-graph

		if [[ -s "${_polap_var_mtcontigs_8mtcontigname}" ]]; then
			_polap_log0_cat "${_polap_var_mtcontigs_8mtcontigname}"
		else
			_polap_log0 "  no such file: ${_polap_var_mtcontigs_8mtcontigname}"
		fi
	done

	# Create a temporary file to store hashes and file paths
	local hash_file="${_polap_var_ga_mtcontigs}/file_hashes.tmp"
	>"$hash_file" # Ensure the file is empty
	for i in "${knum_array[@]}"; do
		KNUM=$i
		source "$script_dir/polap-variables-mtcontigs.sh"
		if [[ -s "${_polap_var_mtcontigs_8mtcontigname}" ]]; then
			md5sum "${_polap_var_mtcontigs_8mtcontigname}" >>"$hash_file"
		fi
	done

	if [[ -s "${hash_file}" ]]; then
		_polap_log2_cat "${hash_file}"
	else
		die "ERROR: no seed contig file of the contig selection types!"
	fi

	# Sort by hash and remove duplicates, keeping only the first occurrence
	# This will give a list of unique files by content
	awk '!seen[$1]++' "$hash_file" >"${_polap_var_ga_mtcontigs}/unique_files_by_content.txt"

	# Display the result
	_polap_log0 "List of unique files by content:"
	_polap_log0 "mt.contig.name files: ${_polap_var_ga_mtcontigs}/unique_files_by_content.txt"

	# Temporary file to store files with fewer than 10 lines
	filtered_file="${_polap_var_ga_mtcontigs}/filtered_files.tmp"
	>"$filtered_file" # Ensure the file is empty

	rm -f "${_polap_var_ga}"/mt.contig.name-*

	# Loop through each unique file
	while IFS=" " read -r hash file; do
		# Compute the new filename based on JNUM and index
		# local new_filename="${_polap_var_ga}/mt.contig.name-${index}"
		# local new_filename_table="${_polap_var_ga}/mtcontig-annotation-table-seed-${index}.txt"

		# Check if file is not empty and line count is less than 10
		if [[ -f "$file" ]]; then
			line_count=$(wc -l <"$file")
			if ((line_count < 30)); then
				# Copy the unique file to the new filename
				# _polap_log3_cmd cp "$file" "$new_filename"
				# third_field=$(echo "$file" | awk -F '/' '{print $4}')
				# cp "${_polap_var_ga_mtcontigs}/${third_field}/8-mtcontig-annotation-table-seed.txt" \
				# 	"${new_filename_table}"
				echo "$line_count $file" >>"$filtered_file"
			fi
		fi
		# Increment the index
		# ((index++))
	done <"${_polap_var_ga_mtcontigs}/unique_files_by_content.txt"

	# Initialize index counter
	local index=${JNUM}

	# Loop through each unique file
	while IFS=" " read -r line_count file; do
		# Compute the new filename based on JNUM and index
		local new_filename="${_polap_var_ga}/mt.contig.name-${index}"
		local new_filename_table="${_polap_var_ga}/mtcontig-annotation-table-seed-${index}.txt"

		# Check if file is not empty and line count is less than 10
		if [[ -f "$file" ]]; then
			line_count=$(wc -l <"$file")
			# Copy the unique file to the new filename
			_polap_log3_cmd cp "$file" "$new_filename"
			third_field=$(echo "$file" | awk -F '/' '{print $4}')
			cp "${_polap_var_ga_mtcontigs}/${third_field}/8-mtcontig-annotation-table-seed.txt" \
				"${new_filename_table}"
		fi
		# Increment the index
		((index++))
	done <"${filtered_file}"

	_polap_log0 "  seed contig name files:"
	ls "${_polap_var_ga}"/mt.contig.name-* >&3

	for file in "${_polap_var_ga}"/mt.contig.name-*; do
		# Extract the <number> part using parameter expansion
		file=$(basename $file)
		local number="${file#mt.contig.name-}"

		_polap_log0 "NEXT: $0 assemble2 -o ${ODIR} -i ${INUM} -j ${number}"
	done
	# if [[ -s "$filtered_file" ]]; then
	# 	# Select the file with the largest number of lines from the filtered files
	# 	local largest_file=$(sort -nr "$filtered_file" | head -n 1 | cut -d ' ' -f 2-)
	#
	# 	# Display the result
	# 	_polap_log0 "    File with unique content and largest number of lines (under 10 lines): ${largest_file}"
	# 	_polap_log3_pipe "cp ${largest_file} ${_polap_var_mtcontigname}"
	# else
	# 	_polap_log0 "ERROR: no seed contig files with less than 10 contigs."
	# fi

	_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Selects contigs using gene density and custom depth range.
#
################################################################################
function _run_polap_seeds-graph() { # select seed contigs
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh"
	source "$script_dir/polap-variables-mtcontigs.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Select contigs based on gene density, a customized depth of coverage, and 
# an analysis of the assembly graph.
#
# The working folder is ${_polap_var_mtcontigs}. Whenever we invoke this menu, 
# the working folder is automatically initialized. We have two depth-range:
# 1. a depth range for the contig preselection.
# 2. another depth range for the graph filtering in later steps.
# The two depth ranges used in the automated process are identical.
#
# TODO:
# 1. document automatic depth-range in step 1
# 2. document preselection in step 2
# 3. document run-polap-r-final-mtcontig.R in step 8
# FIXME: contig version needs to be replaced by edge version.
#
# 1: determine the depth-range either manually or automatically
#   Determining the depth range can be done either by manual adjustment or 
#   through an automated process.
#
# 2: pre-select contigs based on organelle gene annotation
#   Select contigs based on their association with organelle-specific genes
#   for further graph-based filtering.
#
# 3: filter out node sequences of the graph_final.gfa GFA using the depth range
#   Filtering the GFA based on a specified depth range to exclude contigs 
#   that fall outside of this range.
#   3-1: create a no-sequence gfa
#   3-2: extract node sequence part
#   3-3: filter node sequences by the graph-filter depth
#
# 4: prepare the graph for finding connected components before step 5
#   Preparation of the graph is a crucial step in identifying its 
#   connected components.
#   4-1: filter the graph_final.gfa by the depth-range (somewhat complicated)
#   4-2: extract links from the depth-filtered gfa
#   4-3: prepare graphs of the links prior to the step 5
#     links_order: graph index number and edge_number
#
# 5: find connected components by the depth-filtered contigs using DFS algorithm
#   Identify connected components using depth-filtered contigs.
#   links_number: graph
#   links_contig: starting contig for connected components
#   links_seed: graph index numbers of the connected components
#
# 6: convert the connected components to contigs in edge with numbers
#   Convert the depth-filtered connected components into a format that
#   we could use in POLAP organelle-genome assembly.
#   Convert the number-format seed contigs to one in the mt.contig.name format
#   links_seed: graph index numbers of the connected components
#   links_order: graph index number and edge_number
#   links_mtcontig: seed in mt.contig.name format
#
# 7: concatenate the graph depth-filtered edges and annotation-filtered edges
#   Some edges do not belong to the graph because they are preselected but
#   filtered out in the graph-filtering. This happens because we have two
#   kinds of depth ranges: one for preselection and another for graph-filtering.
#   Preselected edges should be in the final short list. However, they could be
#   filtered out in the graph-filtering because we could use a different
#   depth range in the graph-filtering than the preselection depth range.
#   So, we should not have this situation in the automatic seed selection.
#   They are needed to be part of the final mt.contig.name.
#
# 8: rearrange the output files: .table.tsv and filter by max length 1 Mb of 
#   contig seeds
#   This is somewhat weired because we are trying to filter out some seed
#   contigs based on the lengths.
#   It is because the graph-based selection is not perfect.
#   We could filter out too long contigs that are unlikely to be seeds.
#   Graph-filtered seed contigs could be very long because we select contigs
#   that are connected to the preselected seed contigs. We do not consider 
#   seed contig lengths at all.
#   After all, we re-filter the graph-filtered result once again.
#   Here, r-final-mtcontig.R is used.
#   r-final-mtcontig.R also check-marks contigs that are selected as seeds
#   in the annotation depth tabel.
#   report the selected seed contigs
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   --plastid
# Inputs:
#   ${_polap_var_assembly_graph_final_gfa}
#   ${_polap_var_annotation_table}
#   ${_polap_var_wga}/1-mtcontig.depth.stats.txt for manual depth range in --select-contig 1
# Outputs:
#   ${_polap_var_mtcontigs_7mtcontigname}
# Menus:
#   auto -> creating the depth-range values, then making a mt.contig.name
#   manual
#   add -> adding a contig directly to the mt.contig.name
#   remove -> removing a contig directly from the mt.contig.name
#   bandage -> creates a mt.contig.name from selection of contigs in bandage
# View:
#   depth-range <- for step 1
#   preselection <- for step 2
Example: $0 ${_arg_menu[0]} -k 1
Example: $0 ${_arg_menu[0]} -k 2
Example: $0 ${_arg_menu[0]} manual --plastid <- for test (range: 4~9)
Example: $0 ${_arg_menu[0]} manual 1 <- two depth-range files are the same
Example: $0 ${_arg_menu[0]} manual 2 <- two different depth-range files
Example: $0 ${_arg_menu[0]} manual two-value-textfile.txt (for both depth-range)
Example: $0 ${_arg_menu[0]} add|remove
Example: $0 ${_arg_menu[0]} bandage
Example: $0 ${_arg_menu[0]} add new <- for a new mt.contig.name
Example: $0 ${_arg_menu[0]} view -k 1
Example: $0 ${_arg_menu[0]} view table -k 2
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		# for step 1
		if [[ "${_arg_menu[2]}" == "depth-range" ]]; then
			_polap_log0_column "${_polap_var_mtcontigs_depth_range_preselection}"
			_polap_log0_column "${_polap_var_mtcontigs_depth_range_graphfilter}"
		fi
		# for step 2
		if [[ "${_arg_menu[2]}" == "preselection" ]]; then
			_polap_log0_column "${_polap_var_mtcontigs_preselection}"
		fi

		if [[ "${_arg_menu[2]}" == "table" ]]; then
			_polap_log0 "---"
			_polap_log0_cat "${_polap_var_mtcontigs_1_custom_depth_range}"
			local MTCONTIGNAME="${_polap_var_ga}"/mt.contig.name-"${_arg_menu[2]}"
			# cat "${_polap_var_ga}"/mt.contig.name-"${_arg_menu[2]}" >&3
			if [[ -s "${MTCONTIGNAME}" ]]; then
				_polap_log1_cat "${MTCONTIGNAME}"
				_polap_log0 "---"
				if [[ "${_arg_log_stderr}" = "off" ]]; then
					paste -sd',' "${MTCONTIGNAME}" >&3
				else
					paste -sd',' "${MTCONTIGNAME}" >&2
				fi
			fi
			if [[ -s "${_polap_var_mtcontig_table}" ]]; then
				_polap_log0_column "${_polap_var_mtcontig_table}"
				_polap_log0 "---------------------------------"
			else
				_polap_log0 "No such file: ${_polap_var_mtcontig_table}"
			fi

			if [[ -s "${_polap_var_ga_annotation_depth_table_seed_target}" ]]; then
				_polap_log0_column "${_polap_var_ga_annotation_depth_table_seed_target}"
				_polap_log0 "---------------------------------"
			else
				_polap_log0 "No such file: ${_polap_var_ga_annotation_depth_table_seed_target}"
			fi
		else
			# Check if there are no files matching the pattern
			if [ -z "$(ls -A ${_polap_var_ga}/mt.contig.name-*)" ]; then
				_polap_log0 "No seed contigs, yet!"
			else
				if [[ "${_arg_log_stderr}" = "off" ]]; then
					wc -l "${_polap_var_ga}"/mt.contig.name-* >&3
				else
					wc -l "${_polap_var_ga}"/mt.contig.name-* >&2
				fi
			fi
		fi

		return
	fi

	# Add or remove seed contigs
	if [[ "${_arg_menu[1]}" == "add" ]]; then

		if [[ "${_arg_menu[2]}" == "new" ]]; then
			_polap_log0 "you are making a completely new seed contig file: ${MTCONTIGNAME}"
			rm -f "${MTCONTIGNAME}"
		fi

		# Function to append edge_ to the number and write to file
		add_edge() {
			read -p "Enter an edge number (or type 'quit' or 'exit' to stop): " input

			# If the input is not a number, check if it's 'quit' or 'exit'
			if [[ "$input" =~ ^[0-9]+$ ]]; then
				# Prepend "edge_" to the number
				edge_name="edge_$input"

				# Append the result to the file (replace with your file path)
				echo "$edge_name" >>"${MTCONTIGNAME}"

				echo "Appended $edge_name to the file."
			elif [[ "$input" == "quit" || "$input" == "exit" ]]; then
				echo "Exiting and removing duplicates from the file..."

				# Remove duplicates by sorting and keeping only unique entries
				sort -u ${MTCONTIGNAME} -o ${MTCONTIGNAME}
				echo "Duplicates removed. Goodbye!"

				# _polap_seeds_final-mtcontig

				exit

			else
				echo "Invalid input. Please enter a valid number or type 'quit'/'exit' to stop."
			fi
		}

		# Main loop to continuously ask for numbers
		while true; do
			add_edge
		done

		return $RETURN_SUCCESS
	elif [[ "${_arg_menu[1]}" = "bandage" ]]; then
		# Prompt the user for input
		read -p "Enter edges using Bandage (e.g., edge_265, edge_520, edge_425): " edges

		# Convert the input string into new lines (replace ", " with newlines)
		local formatted_edges=$(echo "$edges" | tr ', ' '\n')
		echo "${formatted_edges}" >"${MTCONTIGNAME}"
		# _polap_seeds_final-mtcontig
		return $RETURN_SUCCESS
	fi

	# We initiate the process of selecting seed contigs.
	_polap_log0 "selecting seed contigs using the assembly graph: seed selection type=${KNUM} of ${INUM} (source) -> ${JNUM} (target) ..."
	_polap_log1 "  input1: ${_polap_var_ga_contigger_edges_gfa}"
	_polap_log1 "  input2: ${_polap_var_annotation_table}"

	check_file_existence "${_polap_var_ga_contigger_edges_gfa}"
	check_file_existence "${_polap_var_annotation_table}"

	# step 1. manual custom depth-range or automatic depth-range
	_polap_log1 "  step 1: determining the depth-range either manually or automatically ..."

	# A fresh start at the mtcontigs folder.
	_polap_log2 "  cleaning up (delete and create) the base mtcontigs folder: ${_polap_var_mtcontigs}"
	_polap_log3_cmd rm -rf "${_polap_var_mtcontigs}"
	_polap_log3_cmd mkdir -p "${_polap_var_mtcontigs}"

	_polap_log2 "    output1: ${_polap_var_mtcontigs_depth_range_graphfilter}"
	_polap_log2 "    output2: ${_polap_var_mtcontigs_depth_range_graphfilter}"

	local _is_auto="off"
	if [[ "${_arg_menu[1]}" == "manual" ]]; then
		# manual: menu3 could be a two-value text file.
		_polap_log1 "    create one or two depth-range text files manually ..."
		_polap_log2 "    output1: ${_polap_var_mtcontigs_1_custom_depth_range}"
		_polap_log2 "    output2: ${_polap_var_mtcontigs_2_custom_depth_range}"
		_polap_seeds_create-manual-depth-range "${_arg_menu[2]}"
	else
		_is_auto="on"
		_polap_log1 "  determines the depth range automatically ..."
		_polap_log2 "    input1: ${_polap_var_annotation_table}"
		_polap_log2 "    output1: ${_polap_var_mtcontigs_2_depth_range_by_cdf_copy_number}"

		local _returned_result=""
		_polap_seeds_create-automatic-depth-range "${KNUM}" _returned_result
		if [[ "${_returned_result}" != "depth-range" ]]; then
			return 0
		fi
	fi

	_polap_log3_pipe "cp ${_polap_var_mtcontigs_1_custom_depth_range} \
    ${_polap_var_mtcontigs_depth_range_preselection}"
	_polap_log3_pipe "cp ${_polap_var_mtcontigs_2_custom_depth_range} \
    ${_polap_var_mtcontigs_depth_range_graphfilter}"

	# step 2. preselect contigs for starting contigs later in graph-filtering
	# seed-contig selection.
	_polap_log1 "  step 2: pre-select contigs based on organelle gene annotation"
	_polap_log2 "    input1: ${_polap_var_mtcontigs_depth_range_preselection}"
	_polap_log2 "    input2: ${_polap_var_annotation_table}"
	_polap_log2 "    output: ${_polap_var_mtcontigs_preselection}"
	_polap_seeds_preselect-contigs

	# Continue for the select contig types 3, 4, or 5
	_polap_log1 "  step 3: filtering GFA using the depth range"
	_polap_log2 "    input1: ${_polap_var_assembly_graph_final_gfa}"
	_polap_log2 "    input2: ${_polap_var_mtcontigs_depth_range_graphfilter}"
	_polap_seeds_depthfilter-gfa
	_polap_log2 "    output: ${_polap_var_mtcontigs_gfa_seq_filtered}"

	# Recreate GFA based on filtered edge sequences.
	_polap_log1 "  step 4: preparing the graph for finding connected components"
	_polap_log2 "    input1: ${_polap_var_mtcontigs_gfa_seq_filtered}"
	_polap_log2 "    output1: ${_polap_var_mtcontigs_links_number}"
	_polap_log2 "    output2: ${_polap_var_mtcontigs_links_order}"
	_polap_log2 "    output3: ${_polap_var_mtcontigs_links_contig}"
	_polap_log2 "    output4: ${_polap_var_mtcontigs_links_contig_na}"
	_polap_seeds_prepare-cc
	if [[ ! -s "${_polap_var_mtcontigs_links_tsv}" ]]; then
		_polap_log0 "    no links; so stop here just the pre-selected contigs"
		return 0
	fi

	# Find connected components using Python script
	_polap_log1 "  step 5: finding connected components by the depth-filtered contigs"
	_polap_log2 "    input1: ${_polap_var_mtcontigs_links_number}"
	_polap_log2 "    input2: ${_polap_var_mtcontigs_links_contig}"
	_polap_log2 "    output: ${_polap_var_mtcontigs_links_seed}"
	_polap_log3_pipe "python $script_dir/run-polap-py-find-cc.py \
		${_polap_var_mtcontigs_links_number} \
		${_polap_var_mtcontigs_links_contig} \
		${_polap_var_mtcontigs_links_seed} \
		2>$_polap_output_dest"

	if [[ -s "${_polap_var_mtcontigs_links_seed}" ]]; then
		_polap_log2 "    seed links are found:"
		_polap_log2_head "${_polap_var_mtcontigs_links_seed}"
	fi

	# Choose final mitochondrial contigs
	_polap_log1 "  step 6: converting the depth-filtered contigs in edge with numbers"
	_polap_log2 "    input1: ${_polap_var_mtcontigs_links_seed}"
	_polap_log2 "    input2: ${_polap_var_mtcontigs_links_order}"
	_polap_log2 "    output: ${_polap_var_mtcontigs_links_mtcontig}"
	_polap_log3_pipe "Rscript $script_dir/run-polap-r-cc2mtcontig.R \
		--seed ${_polap_var_mtcontigs_links_seed} \
		--order ${_polap_var_mtcontigs_links_order} \
		--out ${_polap_var_mtcontigs_links_mtcontig} \
		2>$_polap_output_dest"

	if [[ -s "${_polap_var_mtcontigs_links_mtcontig}" ]]; then
		_polap_log3_head "${_polap_var_mtcontigs_links_mtcontig}"
	fi

	_polap_log2 "  step 7: concatenating the graph depth-filtered edges and annotation-filtered NA edges: ${MTCONTIGNAME}"
	_polap_log2 "    input1: ${_polap_var_mtcontigs_links_mtcontig}"
	_polap_log2 "    input2: ${_polap_var_mtcontigs_links_contig_na}"
	_polap_log2 "    output: ${_polap_var_mtcontigs_7mtcontigname}"
	cat "${_polap_var_mtcontigs_links_mtcontig}" "${_polap_var_mtcontigs_links_contig_na}" |
		sort | uniq >"${_polap_var_mtcontigs_7mtcontigname}"

	if [[ -s "${_polap_var_mtcontigs_links_contig_na}" ]]; then
		_polap_log2_head "${_polap_var_mtcontigs_links_contig_na}"
	fi

	_polap_seeds_final-mtcontig "${_is_auto}"

	if [[ -s "${_polap_var_mtcontig_table}" ]]; then
		_polap_log3_column "${_polap_var_mtcontig_table}"
		_polap_log2_pipe "cut -f1 ${_polap_var_mtcontig_table} \
      >${_polap_var_mtcontigs_8mtcontigname}"
	else
		die "ERROR: no final mtcontig: ${_polap_var_mtcontig_table}"
	fi

	if [[ -s "${_polap_var_mtcontigs_8mtcontigname}" ]]; then
		_polap_log1_cat "${_polap_var_mtcontigs_8mtcontigname}"
		_polap_log0 "---"
		if [[ "${_arg_log_stderr}" = "off" ]]; then
			paste -sd',' "${_polap_var_mtcontigs_8mtcontigname}" >&3
		else
			paste -sd',' "${_polap_var_mtcontigs_8mtcontigname}" >&2
		fi
	fi

	local ANUMNEXT=$((INUM + 1))
	_polap_log1 NEXT: $0 map-reads -o "$ODIR" [-i $INUM] [-j $ANUMNEXT]
	_polap_log1 NEXT: $0 assemble2 -o "$ODIR" [-i $INUM] [-j $ANUMNEXT]

	_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Selects contigs for an organelle-genome assembly.
#
# 1. We could select mitochondrial- or plastid-derived contigs using a contig annotation table.
# 2. We determine the range of sequencing depths for those candidate contigs: mean +/- sd \* 3.
# 3. For a given gfa of a genome assembly graph, subset the graph for selecting graph elements in the range.
# 4. Determine connected components in the subset.
# 5. Choose connected components with candidate edges.
#
# We need to read GFA files to manipulate.
# We need to determine connected components.
################################################################################
function _run_polap_seeds-gene() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-common.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Selects contigs using gene density.
# NOTE: a copy of _run_polap_seeds
#
# ------------------------------------------------------------------------------
# Two manual selection approach
# Preselection
# 1: annotation only
# 2: annotation + depth range of depth cumulative distribution
# 3: annotation + depth range of depth mixture distribution
# 4: annotation + custom depth range
# Graph-based selection
# 5: 4 + graph filtered by depth range of depth cumulative distribution
# 6: 4 + graph filtered by depth range of depth mixture distribution
# 7: 4 + graph filtered by custom depth range
#
# Annotation or Gene density 
# Annotation. Depth range. 
# Depth range
# 1. Manual
# 2. Total length ordered by multiplicity
# 3. Multiplicity mixture distribution 
#
# Graph
# Annotation. Depth range. 
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#   -j $JNUM: destination Flye organelle assembly number
#   --select-contig [1|2|3]
# Inputs:
#   ${_polap_var_assembly_graph_final_gfa}
#   ${_polap_var_annotation_table}
#   ${_polap_var_wga}/1-mtcontig.depth.stats.txt for manual depth range in --select-contig 1
# Outputs:
#   ${MTCONTIGNAME}
Example: $0 ${_arg_menu[0]} -i <arg> -j <arg> [auto]
Example: $0 ${_arg_menu[0]} -i <arg> -j <arg> <only|annotation>
Example: $0 ${_arg_menu[0]} -i <arg> -j <arg> <manual|depth>
Example: $0 ${_arg_menu[0]} -i <arg> -j <arg> view
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [[ -s "${_polap_var_mtcontigs_preselection}" ]]; then
			_polap_log0_column "${_polap_var_mtcontigs_preselection}"
		fi
		if [[ -s "${_polap_var_ga_annotation_depth_table_seed_target}" ]]; then
			_polap_log0_column "${_polap_var_ga_annotation_depth_table_seed_target}"
			_polap_log0 "---------------------------------"
		else
			_polap_log0 "No such file: ${_polap_var_ga_annotation_depth_table_seed_target}"
		fi
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	# Two depth-range:
	# 1. depth range for the contig preselection in the step 1.
	# 2. depth range for the graph filtering in later steps.

	_polap_log0 "selecting seed contigs using genes and depths: ${INUM} (source) -> ${JNUM} (target) ..."
	_polap_log1 "  cleaning up (delete and create) the base mtcontigs folder: ${_polap_var_mtcontigs}"
	_polap_log3_cmd rm -rf "${_polap_var_mtcontigs}"
	_polap_log3_cmd mkdir -p "${_polap_var_mtcontigs}"
	check_file_existence "${_polap_var_assembly_graph_final_gfa}"
	check_file_existence "${_polap_var_annotation_table}"
	_polap_log1 "  input1: ${_polap_var_assembly_graph_final_gfa}"
	_polap_log1 "  input2: ${_polap_var_annotation_table}"
	_polap_log1 "  input3: ${_polap_var_mtcontigs_1_custom_depth_range}"
	_polap_log1 "  input4: ${_polap_var_mtcontigs_2_custom_depth_range}"

	# Two depth-range:
	# 1. depth range for the contig preselection in the step 1.
	# 2. depth range for the graph filtering in later steps.

	# Manual custom depth-range or automatic depth-range
	_polap_log1 "  ---"
	_polap_log1 "  step 1: determining the depth-range either manually or automatically ..."
	local _is_auto="off"
	if [[ "${_arg_menu[1]}" == "manual" ]] || [[ "${_arg_menu[1]}" == "depth" ]]; then
		if [[ -s "${_arg_menu[2]}" ]]; then
			_polap_log1 "    input3 (two-value file): ${_arg_menu[2]}"
			_polap_log2_cat "${_arg_menu[2]}"
			numbers=$(grep -Eo '[0-9]+' "${_arg_menu[2]}")
			_polap_log1 "    manual depth range: ${array[0]} ~ ${array[1]}"
			_polap_log3_cmd bash "$script_dir/run-polap-sh-create-depth-file.sh" "${array[0]}" "${array[1]}" "${_polap_var_mtcontigs_1_custom_depth_range}"
			_polap_log2_cat "${_polap_var_mtcontigs_1_custom_depth_range}"
			_polap_log1 "    use this depth-range for both contig preselection and graph filtering"
			_polap_log3_cmd cp "${_polap_var_mtcontigs_1_custom_depth_range}" "${_polap_var_mtcontigs_2_custom_depth_range}"
		else
			# 1. for the depth-range in the contig preselection
			_polap_log0 "  Choose desired values of the depth range for the contig preselection:"
			_polap_log0 "    we will use the depth-range for the graph filtering if you do not have another depth-range file."
			_polap_log0 "    1. First, calculate an estimated average value accordingly."
			_polap_log0 "    2. Choosing a lower bound that is one-third of the average value can be a suitable option."
			_polap_log0 "    3. Selecting an appropriate upper limit can be achieved by setting it to three times the average."
			read -p "  Enter the lower bound of a depth range: " _num1
			read -p "  Enter the upper bound of a depth range: " _num2
			_polap_log3_cmd bash "$script_dir/run-polap-sh-create-depth-file.sh" "${_num1}" "${_num2}" "${_polap_var_mtcontigs_1_custom_depth_range}"
			_polap_log0 "    manual depth range: ${_num1} ~ ${_num2}"
			_polap_log2_cat "${_polap_var_mtcontigs_1_custom_depth_range}"
			_polap_log1 "    use this depth-range for both contig preselection and graph filtering"
			_polap_log3_cmd cp "${_polap_var_mtcontigs_1_custom_depth_range}" "${_polap_var_mtcontigs_2_custom_depth_range}"
		fi
	elif [[ "${_arg_menu[1]}" == "only" ]] || [[ "${_arg_menu[1]}" == "annotation" ]]; then
		_polap_log3_cmd bash "$script_dir/run-polap-sh-create-depth-file.sh" "0" "0" "${_polap_var_mtcontigs_1_custom_depth_range}"
		_polap_log2_cat "${_polap_var_mtcontigs_1_custom_depth_range}"
		_polap_log1 "    use this depth-range for both contig preselection and graph filtering"
		_polap_log3_cmd cp "${_polap_var_mtcontigs_1_custom_depth_range}" "${_polap_var_mtcontigs_2_custom_depth_range}"
	else
		_is_auto="on"
		_polap_log1 "    determines the depth range automatically ..."
		_polap_log2 "      using the length distribution of contigs sorted"
		_polap_log2 "      by copy numbers and organelle gene counts"
		_polap_log2 "      input1: ${_polap_var_annotation_table}"
		_polap_log2 "      output1: ${_polap_var_mtcontigs_2_depth_range_by_cdf_copy_number}"
		local _command1="Rscript $script_dir/run-polap-r-determine-depth-range.R \
				-t ${_polap_var_annotation_table} \
				-o ${_polap_var_mtcontigs_2_depth_range_by_cdf_copy_number}"
		if [[ "${_arg_plastid}" = "on" ]]; then
			_command1+=" \
        --plastid"
		fi
		_command1+=" \
				2>$_polap_output_dest"
		_polap_log3_pipe "${_command1}"

		if [[ -s "${_polap_var_mtcontigs_2_depth_range_by_cdf_copy_number}" ]]; then
			_polap_log3_cat "${_polap_var_mtcontigs_2_depth_range_by_cdf_copy_number}"
			_polap_log1 "    we use one depth-range for contig preselection and graph filtering"
			_polap_log3_cmd cp "${_polap_var_mtcontigs_2_depth_range_by_cdf_copy_number}" "${_polap_var_mtcontigs_1_custom_depth_range}"
			_polap_log3_cmd cp "${_polap_var_mtcontigs_2_depth_range_by_cdf_copy_number}" "${_polap_var_mtcontigs_2_custom_depth_range}"
		else
			die "No automatically generated depth range: no such file: ${_polap_var_mtcontigs_2_depth_range_by_cdf_copy_number}"
		fi
	fi

	# select seed contig preselection

	_polap_log1 "  step 2: select contigs based on organelle gene annotation"

	if [[ -s "${_polap_var_mtcontigs_1_custom_depth_range}" ]]; then
		# Extract numbers from the file
		local numbers=$(awk 'NR==2 {print $1, $2}' "${_polap_var_mtcontigs_1_custom_depth_range}")
		# Store the numbers in variables
		local depth_lower=$(echo $numbers | cut -d ' ' -f 1)
		local depth_upper=$(echo $numbers | cut -d ' ' -f 2)
		_polap_log0 "  depth range for edge contig preselection: $depth_lower ~ $depth_upper"
		_polap_log2_cat "${_polap_var_mtcontigs_1_custom_depth_range}"
	else
		die "ERROR: no such file: ${_polap_var_mtcontigs_1_custom_depth_range}"
		local depth_lower=0
		local depth_upper=0
		_polap_log0 "    custom depth range for starting edge contigs: $depth_lower ~ $depth_upper"
	fi

	_polap_log2 "    minimum gene density for mtDNA: 10 per 1 Mb"
	_polap_log2 "    gene count comparison: MT > PT"
	_polap_log2 "    depth range: $depth_lower ~ $depth_upper"
	_polap_log2 "    input1: ${_polap_var_annotation_table}"
	_polap_log2 "    output: ${_polap_var_mtcontigs_preselection}"
	local _command1="Rscript $script_dir/run-polap-r-preselect-annotation.R  \
		--table ${_polap_var_annotation_table} \
		--out ${_polap_var_mtcontigs_preselection} \
		--depth-range ${depth_lower},${depth_upper} \
		--compare-mt-pt \
		--gene-density 10"
	if [[ "${_arg_plastid}" = "on" ]]; then
		_polap_log2 "    for plastid not mitochondrial DNA"
		_command1+=" \
        --plastid"
	fi
	_command1+=" \
				2>$_polap_output_dest"
	_polap_log3_pipe "${_command1}"

	if [[ -s "${_polap_var_mtcontigs_preselection}" ]]; then
		_polap_log2_cat "${_polap_var_mtcontigs_preselection}"
	else
		_polap_log0 "ERROR: no such file: ${_polap_var_mtcontigs_preselection}"
	fi

	if [[ -s "${_polap_var_mtcontigs_preselection}" ]]; then
		_polap_log0_column "${_polap_var_mtcontigs_preselection}"
		_polap_log2 "  creating the seed contig file: ${MTCONTIGNAME}"
		_polap_log2 "    input1: ${_polap_var_mtcontigs_preselection}"
		_polap_log2 "    output1: ${MTCONTIGNAME}"
		_polap_log3_pipe "cut -f1 ${_polap_var_mtcontigs_preselection} >${MTCONTIGNAME}"
	else
		_polap_log0 "  no such file: ${_polap_var_mtcontigs_preselection}"
	fi

	_polap_log0 "  output1: ${MTCONTIGNAME}"

	_polap_seeds_final-mtcontig

	_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
