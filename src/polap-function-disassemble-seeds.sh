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
# _polap_disassemble_seeds_get-depth-range-of depth-range.txt depth_range
# ${depth_raneg[0]}
# ${depth_raneg[1]}
function _polap_disassemble_seeds_get-depth-range-of {
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
function _polap_disassemble_seeds_create-manual-depth-range {
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
# determines the depth range automatically ..."
# using the length distribution of contigs sorted"
# by copy numbers and organelle gene counts"
################################################################################
_polap_disassemble_seeds_create-automatic-depth-range() {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	local _type=$1
	local -n result_ref=$2
	local _ga_annotation_all="$3"
	local _ga_annotation_cdf_table="$4"
	local _mtcontigs_1_custom_depth_range="$5"
	local _mtcontigs_2_custom_depth_range="$6"

	if [[ "${_arg_plastid}" == "on" ]]; then
		_polap_log2 "  using plastid depth range ..."
		local _command1="Rscript $script_dir/run-polap-r-plastid-determine-depth-range_${_type}.R \
				-t ${_ga_annotation_all} \
				-c ${_ga_annotation_cdf_table} \
				-o ${_mtcontigs_1_custom_depth_range} \
        --plastid \
				2>$_polap_output_dest"
	else
		local _command1="Rscript $script_dir/run-polap-r-determine-depth-range_${_type}.R \
				-t ${_ga_annotation_all} \
				-c ${_ga_annotation_cdf_table} \
				-o ${_mtcontigs_1_custom_depth_range} \
				2>$_polap_output_dest"
	fi
	_polap_log3_pipe "${_command1}"

	_polap_log3_column "${_ga_annotation_cdf_table}"

	if [[ -s "${_mtcontigs_1_custom_depth_range}" ]]; then
		_polap_log1 "    we use one depth-range for contig preselection and graph filtering"
		_polap_log3_pipe "cp ${_mtcontigs_1_custom_depth_range} \
      ${_mtcontigs_2_custom_depth_range}"
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
_polap_disassemble_seeds_preselect-contigs() {

	local _mtcontigs_depth_range_preselection="$1"
	local _ga_annotation_all="$2"
	local _mtcontigs_preselection="$3"

	if [[ -s "${_mtcontigs_depth_range_preselection}" ]]; then
		# Extract numbers from the file
		local numbers=$(awk 'NR==2 {print $1, $2}' "${_mtcontigs_depth_range_preselection}")
		# Store the numbers in variables
		local depth_lower=$(echo $numbers | cut -d ' ' -f 1)
		local depth_upper=$(echo $numbers | cut -d ' ' -f 2)
		_polap_log0 "  depth range for edge contig preselection from ${_mtcontigs_depth_range_preselection}: $depth_lower ~ $depth_upper"
		_polap_log2_column "${_mtcontigs_depth_range_preselection}"
	else
		die "ERROR: no such file: ${_mtcontigs_depth_range_preselection}"
	fi

	_polap_log2 "    minimum gene density for mtDNA: 10 per 1 Mb"
	_polap_log2 "    gene count comparison: MT > PT"
	_polap_log2 "    depth range: $depth_lower ~ $depth_upper"
	_polap_log2 "    input1: ${_mtcontigs_depth_range_preselection}"
	_polap_log2 "    input2: ${_ga_annotation_all}"
	_polap_log2 "    output: ${_mtcontigs_preselection}"
	local _command1="Rscript $script_dir/run-polap-r-preselect-annotation.R  \
		--out ${_mtcontigs_preselection} \
		--table ${_ga_annotation_all} \
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

	if [[ -s "${_mtcontigs_preselection}" ]]; then
		_polap_log2_column "${_mtcontigs_preselection}"
	else
		_polap_log0 "ERROR: no such file: ${_mtcontigs_preselection}"
	fi
}

################################################################################
################################################################################
function _polap_disassemble_seeds_depthfilter-gfa {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	local _contigger=$(dirname "$1")
	local _mtcontigs=$(dirname "$2")
	local _ga_contigger_edges_gfa="$1"
	local _mtcontigs_depth_range_graphfilter="$2"
	local _mtcontigs_gfa_seq_filtered="$3"
	local _mtcontigs_gfa_all="${_mtcontigs}/3-gfa.all.gfa"
	local _mtcontigs_gfa_seq_part="${_mtcontigs}/3-gfa.seq.all.tsv"
	local _mtcontigs_gfa_seq_filtered="${_mtcontigs}/3-gfa.seq.depthfiltered.txt"
	local _mtcontigname="${_mtcontigs}/mtcontigname.txt"

	local depth_values=()
	_polap_disassemble_seeds_get-depth-range-of \
		"${_mtcontigs_depth_range_graphfilter}" \
		depth_values

	if [[ "${depth_values[0]}" -gt 0 ]]; then
		local depth_lower="${depth_values[0]}"
		local depth_upper="${depth_values[1]}"
		_polap_log2 "  depth range for graph filtering: $depth_lower ~ $depth_upper"
	else
		die "ERROR: no depth ranges"
	fi

	_polap_log1 "    step 3-1: creating GFA without sequence data: ${_mtcontigs_gfa_all}"
	_polap_log1 "      input1: ${_ga_contigger_edges_gfa}"
	_polap_log2 "      output: ${_mtcontigs_gfa_all}"
	_polap_log3_pipe "gfatools view \
		-S ${_ga_contigger_edges_gfa} \
		>${_mtcontigs_gfa_all} \
		2>$_polap_output_dest"

	_polap_log1 "    step 3-2: extracting sequence part of GFA: ${_mtcontigs_gfa_seq_part}"
	_polap_log2 "      input1: ${_mtcontigs_gfa_all}"
	_polap_log2 "      output: ${_mtcontigs_gfa_seq_part}"
	_polap_log3_pipe "grep ^S ${_mtcontigs_gfa_all} >${_mtcontigs_gfa_seq_part}"

	# Filter edges in GFA using depths.
	_polap_log1 "    step 3-3: filtering GFA sequence part using depth range"
	_polap_log2 "      input1: ${_mtcontigs_gfa_seq_part}"
	_polap_log2 "      input2: ${_mtcontigs_depth_range_graphfilter}"
	_polap_log2 "        depth range: $depth_lower ~ $depth_upper"
	_polap_log2 "      output: ${_mtcontigs_gfa_seq_filtered}"
	if [[ -s "${_mtcontigs_depth_range_graphfilter}" ]]; then
		_polap_log3_pipe "Rscript $script_dir/run-polap-r-depthfilter-gfa.R \
			--gfa ${_mtcontigs_gfa_seq_part} \
			--depth ${_mtcontigs_depth_range_graphfilter} \
			--out ${_mtcontigs_gfa_seq_filtered} \
			2>$_polap_output_dest"
	else
		_polap_log0 "    no such file: ${_mtcontigs_depth_range_graphfilter}"
		_polap_log0 "    ${_mtcontigname} -> empty or no seed contigs"
		>"${_mtcontigname}"
		_polap_log2_cat "${_mtcontigname}"
		_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ -s "${_mtcontigs_gfa_seq_filtered}" ]]; then
		_polap_log3_head "${_mtcontigs_gfa_seq_filtered}"
	fi
}

################################################################################
# Prepare data for connected components
################################################################################
_polap_disassemble_seeds_prepare-cc() {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	local _mtcontigs=$(dirname "$1")
	local _mtcontigs_preselection="$1"
	local _mtcontigs_gfa_seq_filtered="$2"
	local _mtcontigs_links_number="$3"
	local _mtcontigs_links_order="$4"
	local _mtcontigs_links_contig="$5"
	local _mtcontigs_links_contig_na="$6"
	local _mtcontigname="$7"

	local _mtcontigs_gfa_seq_filtered_edge="${_mtcontigs}/4-gfa.seq.depthfiltered.edge.txt"
	local _mtcontigs_gfa_depthfiltered_gfa="${_mtcontigs}/4-gfa.depthfiltered.gfa"
	local _mtcontigs_links="${_mtcontigs}/4-gfa.links"
	local _mtcontigs_links_tsv="${_mtcontigs}/4-gfa.links.tsv"

	_polap_log2 "    step 4-1: subsetting GFA using the depth-filtered GFA sequence part: ${_mtcontigs_gfa_depthfiltered_gfa}"
	_polap_log2 "      input1: ${_mtcontigs_gfa_seq_filtered}"
	_polap_log2 "      output1: ${_mtcontigs_gfa_seq_filtered_edge}"
	_polap_log2 "      output2: ${_mtcontigs_gfa_depthfiltered_gfa}"

	_polap_log3_pipe "cut -f1 \
    ${_mtcontigs_gfa_seq_filtered} \
    >${_mtcontigs_gfa_seq_filtered_edge}"

	_polap_log3_pipe "gfatools view -S \
		-l @${_mtcontigs_gfa_seq_filtered_edge} \
		${_ga_contigger_edges_gfa} \
		2>$_polap_output_dest \
		>${_mtcontigs_gfa_depthfiltered_gfa}"

	if [[ -s "${_mtcontigs_gfa_depthfiltered_gfa}" ]]; then
		_polap_log2_head "${_mtcontigs_gfa_depthfiltered_gfa}"
	fi

	# Prepare links for finding connected components.
	_polap_log2 "    step 4-2: preparing links for finding connected components"
	_polap_log2 "      input: ${_mtcontigs_gfa_depthfiltered_gfa}"
	_polap_log2 "      output: ${_mtcontigs_links_tsv}"
	_polap_log3_pipe "grep ^L ${_mtcontigs_gfa_depthfiltered_gfa} |\
    cut -f2,4 >${_mtcontigs_links_tsv}"

	# No links case
	if [[ -s "${_mtcontigs_links_tsv}" ]]; then
		_polap_log2_head "${_mtcontigs_links_tsv}"
	else
		# no links -> just current selection
		_polap_log3_pipe "cut -f2 ${_mtcontigs_gfa_depthfiltered_gfa} |
			sort | uniq >${_mtcontigname}"
		_polap_log0 "  output1: ${_mtcontigname}"
		_polap_log2_cat "${_mtcontigname}"

		# if [[ "${_arg_plastid}" = "off" ]]; then
		# 	_polap_disassemble_seeds_final-mtcontig \
		# 		"${_ga_annotation_all}" \
		# 		"${_ga_annotation_depth_table}" \
		# 		"${_mtcontigs_7mtcontigname}" \
		# 		"${_mtcontig_table}" \
		# 		"${_mtcontigs_annotation_table_seed}"
		# else
		# 	_polap_disassemble_seeds_final-mtcontig \
		# 		"${_ga_annotation_all}" \
		# 		"${_ga_pt_annotation_depth_table}" \
		# 		"${_mtcontigs_7mtcontigname}" \
		# 		"${_mtcontig_table}" \
		# 		"${_mtcontigs_annotation_table_seed}"
		# fi

		_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	local _preselection="${_mtcontigs_preselection}"

	# Run R script to analyze GFA links
	_polap_log2 "    step 4-3: preparing for finding connected components"
	_polap_log2 "      input1: ${_preselection}"
	_polap_log2 "      input2: ${_mtcontigs_links_tsv}"
	_polap_log2 "      output-base: ${_mtcontigs_links}"
	_polap_log2 "      output1: ${_mtcontigs_links_number}"
	_polap_log2 "      output2: ${_mtcontigs_links_order}"
	_polap_log2 "      output3: ${_mtcontigs_links_contig}"
	_polap_log2 "      output4: ${_mtcontigs_links_contig_na}"
	_polap_log3_pipe "Rscript $script_dir/run-polap-r-prepare-cc.R \
		--edge ${_preselection} \
		--gfa ${_mtcontigs_links_tsv} \
		--out ${_mtcontigs_links} \
		--unit-annotation-edge \
		2>$_polap_output_dest"

	if [[ -s "${_mtcontigs_links_number}" ]]; then
		_polap_log2_head "${_mtcontigs_links_number}"
	fi
	if [[ -s "${_mtcontigs_links_order}" ]]; then
		_polap_log2_head "${_mtcontigs_links_order}"
	fi
	if [[ -s "${_mtcontigs_links_contig}" ]]; then
		_polap_log2_head "${_mtcontigs_links_contig}"
	fi
	if [[ -s "${_mtcontigs_links_contig_na}" ]]; then
		_polap_log2_head "${_mtcontigs_links_contig_na}"
	else
		_polap_log2 "No NA edge: ${_mtcontigs_links_contig_na}"
	fi
}

################################################################################
# Finalize the mt.contig.name.
################################################################################
# TODO
_polap_disassemble_seeds_final-mtcontig() {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	local _ga_annotation_all="$1"
	local _ga_annotation_depth_table="$2"
	local _mtcontigs_7mtcontigname="$3"
	local _mtcontig_table="$4"
	local _mtcontigs_annotation_table_seed="$5"

	# local _is_auto=$1 not using it anymore? for --length 1000000 option below

	_polap_log1 "  step 8-1: filtering by max length 1 Mb of contig seeds"
	_polap_log1 "    1 Mb length filtering for auto"
	_polap_log2 "    input1: ${_ga_annotation_all}"
	_polap_log2 "    input2: ${_mtcontigs_7mtcontigname}"
	_polap_log2 "    output1: ${_mtcontig_table}"
	local _command1="Rscript $script_dir/run-polap-r-final-filter-mtcontig.R \
		-t ${_ga_annotation_all} \
		-m ${_mtcontigs_7mtcontigname} \
    -o ${_mtcontig_table} \
    --length 1000000 \
		2>$_polap_output_dest"
	_polap_log3_pipe "${_command1}"

	_polap_disassemble_seeds_final-seeds-mtcontig \
		"${_ga_annotation_all}" \
		"${_ga_annotation_depth_table}" \
		"${_mtcontigs_7mtcontigname}" \
		"${_mtcontigs_annotation_table_seed}"
}

# TODO
_polap_disassemble_seeds_final-seeds-mtcontig() {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	local _ga_annotation_all="$1"
	local _ga_annotation_depth_table="$2"
	local _mt_contig_name="$3"
	local _annotation_table_seed_target="$4"

	_polap_log1 "  annotation table with contig seed marks"
	_polap_log2 "    input1: ${_ga_annotation_all}"
	_polap_log2 "    input2: ${_ga_annotation_depth_table}"
	_polap_log2 "    input3: ${_mt_contig_name}"
	_polap_log2 "    output1: ${_annotation_table_seed_target}"
	local _command1="Rscript $script_dir/run-polap-r-final-seed-mtcontig.R \
		-t ${_ga_annotation_all} \
		-a ${_ga_annotation_depth_table} \
		-m ${_mt_contig_name} \
    -o ${_annotation_table_seed_target} \
		2>$_polap_output_dest"
	_polap_log3_pipe "${_command1}"

}

function _polap_disassemble_seeds_report-mtcontig {
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

function _run_polap_choose-seed { # select seed contigs
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
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   -k ${_arg_knum}: destination Flye organelle assembly number
#   --plastid
# Inputs:
#   ${_polap_var_mtcontigs_7mtcontigname}
# Outputs:
#   ${_polap_var_mtcontigname}
Example: $(basename "$0") ${_arg_menu[0]} -i 1 -j 2 -k 1
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
	_polap_log0 "choose seed contigs using the assembly graph: ${_arg_inum} (source) -> ${_arg_jnum} (target) ..."

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
# Prepare seed contigs using a Flye genome assembly.
################################################################################
polap_disassemble-seeds() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	# source "$script_dir/polap-variables-common.sh"
	# source "$script_dir/polap-variables-mtcontigs.sh"

	local _contigger=$(dirname "$1")
	local _ga=$(dirname "${_contigger}")
	local _ann="${_ga}/50-annotation"
	local _ga_contigger_edges_gfa="$1"
	local _ga_annotation_all="$2"
	local _mtcontigname="$3"

	_arg_plastid="on"
	local _ga_mtcontigs="${_ga}/51-mtcontigs"
	# local _mtcontigname="${_ga}/mt.contig.name"

	# We initiate the process of selecting seed contigs.
	_polap_log0 "select seed contigs using the assembly graph with multiple selection methods"
	_polap_log1 "  input1: ${_ga_contigger_edges_gfa}"
	_polap_log1 "  input2: ${_ga_annotation_all}"
	_polap_log1 "  output: ${_mtcontigname}"

	rm -f "${_ga}"/mt.contig.name-*

	# Create the array based on the value of _arg_plastid
	if [[ "$_arg_plastid" == "on" ]]; then
		# local knum_array=(1 6)
		local knum_array=(2)
	else
		local knum_array=({1..6}) # This creates an array of 1 through 6
	fi

	# 1. file_hashed.tmp
	# 2. unique_files_by_content.txt
	# 3. filtered_files.tmp

	_polap_log3_cmd mkdir -p "${_ga_mtcontigs}"

	# Create a temporary file to store hashes and file paths
	local hash_file="${_ga_mtcontigs}/file_hashes.tmp"
	>"$hash_file" # Ensure the file is empty
	local i="0"
	for i in "${knum_array[@]}"; do
		_knum=$i
		source "$script_dir/polap-variables-mtcontigs.sh"
		polap_disassemble-seeds-graph \
			"${_ga_contigger_edges_gfa}" \
			"${_ga_annotation_all}" \
			"${_knum}"

		local _mtcontigs="${_ga}/51-mtcontigs/${_knum}"
		local _mtcontigs_8mtcontigname="${_mtcontigs}/8-mt.contig.name.txt"
		if [[ -s "${_mtcontigs_8mtcontigname}" ]]; then
			md5sum "${_mtcontigs_8mtcontigname}" >>"$hash_file"
		fi
	done

	if [[ -s "${hash_file}" ]]; then
		_polap_log3_cat "${hash_file}"
	else
		return "${_POLAP_ERR_NO_SEEDS}"
	fi

	# This will give a list of unique files by content
	awk '!seen[$1]++' "$hash_file" >"${_ga_mtcontigs}/unique_files_by_content.txt"

	# Temporary file to store files with fewer than 10 lines
	filtered_file="${_ga_mtcontigs}/filtered_files.tmp"
	>"$filtered_file" # Ensure the file is empty
	while IFS=" " read -r hash file; do
		# Check if file is not empty and line count is less than _arg_max_seeds
		if [[ -f "$file" ]]; then
			local line_count=$(wc -l <"$file")
			if ((line_count < _arg_max_seeds)); then
				echo "$line_count $file" >>"$filtered_file"
			fi
		fi
	done <"${_ga_mtcontigs}/unique_files_by_content.txt"

	local n=$(wc -l <"$filtered_file")
	if ((n == 1)); then
		IFS=" " read -r line_count file <"$filtered_file"
		_polap_log0 "we have a single seed: $file ($line_count)"
		_polap_log3_cmd cp "$file" "${_mtcontigname}"
	elif ((n == 0)); then
		_polap_log0 "we have no seeds."
	else
		_polap_log0 "we have multiple seeds."
	fi

	return 0

	# NOTE: we do not use the following because plastid genome assemblies tend to have
	# a simple seed conformation.

	# Initialize index counter
	local index=0

	# Loop through each unique file
	while IFS=" " read -r line_count file; do
		# Compute the new filename based on _arg_jnum and index
		local new_filename="${_ga}/mt.contig.name-${index}"
		local new_filename_table="${_ga}/contig-annotation-table-seed-${index}.txt"

		# Check if file is not empty and line count is less than 10
		if [[ -f "$file" ]]; then
			line_count=$(wc -l <"$file")
			# Copy the unique file to the new filename
			_polap_log3_cmd cp "$file" "$new_filename"
			third_field=$(echo "$file" | awk -F '/' '{print $4}')
			cp "${_ga_mtcontigs}/${third_field}/8-mtcontig-annotation-table-seed.txt" \
				"${new_filename_table}"
		fi
		# Increment the index
		((index++))
	done <"${filtered_file}"

	_polap_log0 "seed contig name files:"
	ls "${_ga}"/mt.contig.name-* >&3

	for _mtcontigname in "${_ga}"/mt.contig.name-*; do
		# Extract the <number> part using parameter expansion

		local file=$(basename ${_mtcontigname})
		local number="${file#mt.contig.name-}"
		local _ga_annotation_depth_table_seed_target="${_ga}/contig-annotation-depth-table-seed-${number}.txt"

		_polap_disassemble_seeds_final-seeds-mtcontig \
			"${_mtcontigname}" \
			"${_ga_annotation_depth_table_seed_target}"

		_polap_log0 "NEXT: $(basename "$0") assemble2 -o ${_arg_outdir} -i ${_arg_inum} -j ${number}"
	done

	_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Select contigs using gene density and custom depth range.
################################################################################
polap_disassemble-seeds-graph() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	# source "$script_dir/polap-variables-common.sh"
	# source "$script_dir/polap-variables-mtcontigs.sh"

	local _contigger=$(dirname "$1")
	local _ga=$(dirname "${_contigger}")
	local _ann="${_ga}/50-annotation"
	local _ga_contigger_edges_gfa="$1"
	local _ga_annotation_all="$2"
	local _knum="$3"

	# from polap-variables-mtcontigs.sh
	local _mtcontigs="${_ga}/51-mtcontigs/${_knum}"
	local _mtcontigs_depth_range_preselection="${_mtcontigs}/depth.range.preselection.txt"
	local _mtcontigs_depth_range_graphfilter="${_mtcontigs}/depth.range.graphfilter.txt"
	local _mtcontigs_1_custom_depth_range="${_mtcontigs}/custom.depth.range-1.txt"
	local _mtcontigs_2_custom_depth_range="${_mtcontigs}/custom.depth.range-2.txt"
	local _mtcontigs_preselection="${_mtcontigs}/1-preselection.by.gene.density.txt"
	local _mtcontigs_gfa_all="${_mtcontigs}/3-gfa.all.gfa"
	local _mtcontigs_gfa_seq_filtered="${_mtcontigs}/3-gfa.seq.depthfiltered.txt"
	local _mtcontigs_gfa_seq_part="${_mtcontigs}/3-gfa.seq.all.tsv"
	local _mtcontigs_gfa_seq_filtered_edge="${_mtcontigs}/4-gfa.seq.depthfiltered.edge.txt"
	local _mtcontigs_gfa_depthfiltered_gfa="${_mtcontigs}/4-gfa.depthfiltered.gfa"

	local _mtcontigs_links="${_mtcontigs}/4-gfa.links"
	local _mtcontigs_links_tsv="${_mtcontigs}/4-gfa.links.tsv"
	local _mtcontigs_links_contig_na="${_mtcontigs}/4-gfa.links.contig.na.txt"
	local _mtcontigs_links_contig="${_mtcontigs}/4-gfa.links.contig.txt"
	local _mtcontigs_links_number="${_mtcontigs}/4-gfa.links.number.txt"
	local _mtcontigs_links_order="${_mtcontigs}/4-gfa.links.order.txt"
	local _mtcontigs_links_seed="${_mtcontigs}/5-gfa.links.seed.txt"
	local _mtcontigs_links_mtcontig="${_mtcontigs}/6-gfa.links.mtcontig.txt"
	local _mtcontigs_7mtcontigname="${_mtcontigs}/7-mt.contig.name.txt"
	local _mtcontigs_8mtcontigname="${_mtcontigs}/8-mt.contig.name.txt"
	local _mtcontig_table="${_mtcontigs}/8-mtcontig.table.tsv"
	local _mtcontigs_annotation_table_seed="${_mtcontigs}/8-mtcontig-annotation-table-seed.txt"
	# delete these later
	local _mtcontigs_2_depth_range_by_cdf_copy_number="${_mtcontigs}/2-depth.range.by.cdf.copy.number.txt"
	local _mtcontig_table="${_mtcontigs}/8-mtcontig.table.tsv"
	# delete these later
	local _mtcontigs_2_depth_range_by_cdf_copy_number="${_mtcontigs}/2-depth.range.by.cdf.copy.number.txt"

	# from polap-variables-common.sh
	local _mtcontigname="${_ga}/mt.contig.name"
	local _ga_annotation_depth_table="${_ga}/contig-annotation-depth-table.txt"
	local _ga_annotation_cdf_table="${_ga}/contig-annotation-cdf-table.txt"
	local _ga_pt_annotation_depth_table="${_ga}/pt-contig-annotation-depth-table.txt"

	# We initiate the process of selecting seed contigs.
	_polap_log0 "select seed contigs using the assembly graph"
	_polap_log1 "  input1: ${_ga_contigger_edges_gfa}"
	_polap_log1 "  input2: ${_ga_annotation_all}"

	check_file_existence "${_ga_contigger_edges_gfa}"
	check_file_existence "${_ga_annotation_all}"

	# step 1. manual custom depth-range or automatic depth-range
	_polap_log1 "  step 1: determining the depth-range either manually or automatically ..."

	# A fresh start at the mtcontigs folder.
	_polap_log2 "  cleaning up (delete and create) the base mtcontigs folder: ${_mtcontigs}"
	_polap_log3_cmd rm -rf "${_mtcontigs}"
	_polap_log3_cmd mkdir -p "${_mtcontigs}"

	_polap_log2 "    output1: ${_mtcontigs_depth_range_graphfilter}"

	local _is_auto="on"
	_polap_log1 "  determines the depth range automatically ..."
	_polap_log2 "    input1: ${_ga_annotation_all}"
	_polap_log2 "    output1: ${_mtcontigs_2_depth_range_by_cdf_copy_number}"

	local _returned_result=""
	_polap_disassemble_seeds_create-automatic-depth-range \
		"${_knum}" \
		_returned_result \
		"${_ga_annotation_all}" \
		"${_ga_annotation_cdf_table}" \
		"${_mtcontigs_1_custom_depth_range}" \
		"${_mtcontigs_2_custom_depth_range}"
	if [[ "${_returned_result}" != "depth-range" ]]; then
		return 0
	fi

	_polap_log3_pipe "cp ${_mtcontigs_1_custom_depth_range} \
    ${_mtcontigs_depth_range_preselection}"
	_polap_log3_pipe "cp ${_mtcontigs_2_custom_depth_range} \
    ${_mtcontigs_depth_range_graphfilter}"

	# step 2. preselect contigs for starting contigs later in graph-filtering
	# seed-contig selection.
	_polap_log1 "  step 2: pre-select contigs based on organelle gene annotation"
	_polap_log2 "    input1: ${_mtcontigs_depth_range_preselection}"
	_polap_log2 "    input2: ${_ga_annotation_all}"
	_polap_log2 "    output: ${_mtcontigs_preselection}"
	_polap_disassemble_seeds_preselect-contigs \
		"${_mtcontigs_depth_range_preselection}" \
		"${_ga_annotation_all}" \
		"${_mtcontigs_preselection}"

	# Continue for the select contig types 3, 4, or 5
	_polap_log1 "  step 3: filtering GFA using the depth range"
	_polap_log2 "    input1: ${_ga_contigger_edges_gfa}"
	_polap_log2 "    input2: ${_mtcontigs_depth_range_graphfilter}"
	_polap_log2 "    output: ${_mtcontigs_gfa_seq_filtered}"
	_polap_disassemble_seeds_depthfilter-gfa \
		"${_ga_contigger_edges_gfa}" \
		"${_mtcontigs_depth_range_graphfilter}" \
		"${_mtcontigs_gfa_seq_filtered}"

	# Recreate GFA based on filtered edge sequences.
	_polap_log1 "  step 4: preparing the graph for finding connected components"
	_polap_log2 "    input1: ${_mtcontigs_preselection} new in disassemble-seeds"
	_polap_log2 "    input2: ${_mtcontigs_gfa_seq_filtered}"
	_polap_log2 "    output1: ${_mtcontigs_links_number}"
	_polap_log2 "    output2: ${_mtcontigs_links_order}"
	_polap_log2 "    output3: ${_mtcontigs_links_contig}"
	_polap_log2 "    output4: ${_mtcontigs_links_contig_na}"
	_polap_log2 "    output5: ${_mtcontigname}"
	_polap_disassemble_seeds_prepare-cc \
		"${_mtcontigs_preselection}" \
		"${_mtcontigs_gfa_seq_filtered}" \
		"${_mtcontigs_links_number}" \
		"${_mtcontigs_links_order}" \
		"${_mtcontigs_links_contig}" \
		"${_mtcontigs_links_contig_na}" \
		"${_mtcontigname}"
	if [[ ! -s "${_mtcontigs_links_tsv}" ]]; then
		_polap_log0 "    no links; so stop here just the pre-selected contigs"
		return 0
	fi

	# Find connected components using Python script
	_polap_log1 "  step 5: finding connected components by the depth-filtered contigs"
	_polap_log2 "    input1: ${_mtcontigs_links_number}"
	_polap_log2 "    input2: ${_mtcontigs_links_contig}"
	_polap_log2 "    output: ${_mtcontigs_links_seed}"
	_polap_log3_pipe "python $script_dir/run-polap-py-find-cc.py \
		${_mtcontigs_links_number} \
		${_mtcontigs_links_contig} \
		${_mtcontigs_links_seed} \
		2>$_polap_output_dest"

	if [[ -s "${_mtcontigs_links_seed}" ]]; then
		_polap_log2 "    seed links are found:"
		_polap_log2_head "${_mtcontigs_links_seed}"
	fi

	# Choose final mitochondrial contigs
	_polap_log1 "  step 6: converting the depth-filtered contigs in edge with numbers"
	_polap_log2 "    input1: ${_mtcontigs_links_seed}"
	_polap_log2 "    input2: ${_mtcontigs_links_order}"
	_polap_log2 "    output: ${_mtcontigs_links_mtcontig}"
	_polap_log3_pipe "Rscript $script_dir/run-polap-r-cc2mtcontig.R \
		--seed ${_mtcontigs_links_seed} \
		--order ${_mtcontigs_links_order} \
		--out ${_mtcontigs_links_mtcontig} \
		2>$_polap_output_dest"

	if [[ -s "${_mtcontigs_links_mtcontig}" ]]; then
		_polap_log3_head "${_mtcontigs_links_mtcontig}"
	fi

	_polap_log2 "  step 7: concatenating the graph depth-filtered edges and annotation-filtered NA edges: ${_mtcontigname}"
	_polap_log2 "    input1: ${_mtcontigs_links_mtcontig}"
	_polap_log2 "    input2: ${_mtcontigs_links_contig_na}"
	_polap_log2 "    output: ${_mtcontigs_7mtcontigname}"
	_polap_log2 "    output: ${_mtcontigs_8mtcontigname}"
	cat "${_mtcontigs_links_mtcontig}" "${_mtcontigs_links_contig_na}" |
		sort | uniq >"${_mtcontigs_7mtcontigname}"

	if [[ -s "${_mtcontigs_links_contig_na}" ]]; then
		_polap_log2_head "${_mtcontigs_links_contig_na}"
	fi

	if [[ "${_arg_plastid}" = "off" ]]; then
		_polap_disassemble_seeds_final-mtcontig \
			"${_ga_annotation_all}" \
			"${_ga_annotation_depth_table}" \
			"${_mtcontigs_7mtcontigname}" \
			"${_mtcontig_table}" \
			"${_mtcontigs_annotation_table_seed}"
	else
		_polap_disassemble_seeds_final-mtcontig \
			"${_ga_annotation_all}" \
			"${_ga_pt_annotation_depth_table}" \
			"${_mtcontigs_7mtcontigname}" \
			"${_mtcontig_table}" \
			"${_mtcontigs_annotation_table_seed}"
	fi

	if [[ -s "${_mtcontig_table}" ]]; then
		_polap_log3_column "${_mtcontig_table}"
		_polap_log2_pipe "cut -f1 ${_mtcontig_table} \
      >${_mtcontigs_8mtcontigname}"
	else
		die "ERROR: no final mtcontig: ${_mtcontig_table}"
	fi

	if [[ -s "${_mtcontigs_8mtcontigname}" ]]; then
		_polap_log1_cat "${_mtcontigs_8mtcontigname}"
		_polap_log0 "---"
		if [[ "${_arg_log_stderr}" = "off" ]]; then
			paste -sd',' "${_mtcontigs_8mtcontigname}" >&3
		else
			paste -sd',' "${_mtcontigs_8mtcontigname}" >&2
		fi
	fi

	# Disable debugging if previously enabled
	_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
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
function _run_polap_disassemble_seeds-gene {
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
# NOTE: a copy of _run_polap_disassemble_seeds
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
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   --select-contig [1|2|3]
# Inputs:
#   ${_polap_var_ga_contigger_edges_gfa}
#   ${_polap_var_ga_annotation_all}
#   ${_polap_var_wga}/1-mtcontig.depth.stats.txt for manual depth range in --select-contig 1
# Outputs:
#   ${_polap_var_mtcontigname}
Example: $(basename "$0") ${_arg_menu[0]} -i <arg> -j <arg> [auto]
Example: $(basename "$0") ${_arg_menu[0]} -i <arg> -j <arg> <only|annotation>
Example: $(basename "$0") ${_arg_menu[0]} -i <arg> -j <arg> <manual|depth>
Example: $(basename "$0") ${_arg_menu[0]} -i <arg> -j <arg> view
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

	_polap_log0 "selecting seed contigs using genes and depths: ${_arg_inum} (source) -> ${_arg_jnum} (target) ..."
	_polap_log1 "  cleaning up (delete and create) the base mtcontigs folder: ${_polap_var_mtcontigs}"
	_polap_log3_cmd rm -rf "${_polap_var_mtcontigs}"
	_polap_log3_cmd mkdir -p "${_polap_var_mtcontigs}"
	check_file_existence "${_polap_var_ga_contigger_edges_gfa}"
	check_file_existence "${_polap_var_ga_annotation_all}"
	_polap_log1 "  input1: ${_polap_var_ga_contigger_edges_gfa}"
	_polap_log1 "  input2: ${_polap_var_ga_annotation_all}"
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
		_polap_log2 "      input1: ${_polap_var_ga_annotation_all}"
		_polap_log2 "      output1: ${_polap_var_mtcontigs_2_depth_range_by_cdf_copy_number}"
		local _command1="Rscript $script_dir/run-polap-r-determine-depth-range.R \
				-t ${_polap_var_ga_annotation_all} \
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
	_polap_log2 "    input1: ${_polap_var_ga_annotation_all}"
	_polap_log2 "    output: ${_polap_var_mtcontigs_preselection}"
	local _command1="Rscript $script_dir/run-polap-r-preselect-annotation.R  \
		--table ${_polap_var_ga_annotation_all} \
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
		_polap_log2 "  creating the seed contig file: ${_polap_var_mtcontigname}"
		_polap_log2 "    input1: ${_polap_var_mtcontigs_preselection}"
		_polap_log2 "    output1: ${_polap_var_mtcontigname}"
		_polap_log3_pipe "cut -f1 ${_polap_var_mtcontigs_preselection} >${_polap_var_mtcontigname}"
	else
		_polap_log0 "  no such file: ${_polap_var_mtcontigs_preselection}"
	fi

	_polap_log0 "  output1: ${_polap_var_mtcontigname}"

	_polap_disassemble_seeds_final-mtcontig

	_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
