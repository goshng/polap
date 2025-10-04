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
# Functions for selecting seed contigs that are used for filtering reads for
# organelle-genome assembly.
# I think that manual contig selection is still better than automatic for mtDNA.
# This is under development. Most of the code are used for ptDNA seed selection.
# This may work for some, but I do not recommend its use because it could
# only select seeds that could have been selected by any user. It works for
# some datasets, but it needs more testing.
#
# function _polap_seeds_get-depth-range-of {
# function _polap_seeds_create-manual-depth-range {
# function _polap_seeds_create-automatic-depth-range {
# function _polap_seeds_preselect-contigs {
# function _polap_seeds_depthfilter-gfa {
# function _polap_seeds_prepare-cc {
# function _polap_seeds_final-mtcontig {
# function _polap_seeds_final-seeds-mtcontig {
# function _polap_seeds_report-mtcontig {
# function _run_polap_choose-seed { # select seed contigs
# function _run_polap_seeds { # select seed contigs
# function _run_polap_seeds-graph { # select seed contigs
# function _run_polap_seeds-gene {
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

function _run_polap_choose-seed { # select seed contigs
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"
	source "${_POLAPLIB_DIR}/polap-variables-mtcontigs.sh"

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
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
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
function _run_polap_seeds-gene {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

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
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
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
			_polap_log3_cmd bash "${_POLAPLIB_DIR}/polap-bash-create-depth-file.sh" "${array[0]}" "${array[1]}" "${_polap_var_mtcontigs_1_custom_depth_range}"
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
			_polap_log3_cmd bash "${_POLAPLIB_DIR}/polap-bash-create-depth-file.sh" "${_num1}" "${_num2}" "${_polap_var_mtcontigs_1_custom_depth_range}"
			_polap_log0 "    manual depth range: ${_num1} ~ ${_num2}"
			_polap_log2_cat "${_polap_var_mtcontigs_1_custom_depth_range}"
			_polap_log1 "    use this depth-range for both contig preselection and graph filtering"
			_polap_log3_cmd cp "${_polap_var_mtcontigs_1_custom_depth_range}" "${_polap_var_mtcontigs_2_custom_depth_range}"
		fi
	elif [[ "${_arg_menu[1]}" == "only" ]] || [[ "${_arg_menu[1]}" == "annotation" ]]; then
		_polap_log3_cmd bash "${_POLAPLIB_DIR}/polap-bash-create-depth-file.sh" "0" "0" "${_polap_var_mtcontigs_1_custom_depth_range}"
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
		local _command1="Rscript ${_POLAPLIB_DIR}/polap-r-determine-depth-range.R \
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
	local _command1="Rscript ${_POLAPLIB_DIR}/polap-r-preselect-annotation.R  \
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

	_polap_seeds_final-mtcontig

	_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
