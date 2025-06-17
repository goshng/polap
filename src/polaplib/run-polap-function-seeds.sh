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

# depth_range=()
# _polap_seeds_get-depth-range-of depth-range.txt depth_range
# ${depth_raneg[0]}
# ${depth_raneg[1]}
function _polap_seeds_get-depth-range-of {
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
function _polap_seeds_create-manual-depth-range {
	local _two_value_textfile="${_arg_menu[2]}"

	if [[ -s "${_two_value_textfile}" ]]; then
		_polap_log1 "    input3 (two-value file): ${_two_value_textfile}"
		_polap_log2_cat "${_two_value_textfile}"
		readarray -t _numbers < <(grep -Eo '[0-9]+' "${_two_value_textfile}")
		_polap_log1 "    manual depth range: ${_numbers[0]}"
		_polap_log1 "    manual depth range: ${_numbers[0]} ~ ${_numbers[1]}"
		_polap_log3_cmd bash "${_POLAPLIB_DIR}/run-polap-sh-create-depth-file.sh" "${_numbers[0]}" "${_numbers[1]}" "${_polap_var_mtcontigs_1_custom_depth_range}"
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
		_polap_log3_cmd bash "${_POLAPLIB_DIR}/run-polap-sh-create-depth-file.sh" "${_num1}" "${_num2}" "${_polap_var_mtcontigs_1_custom_depth_range}"
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
			_polap_log3_cmd bash "${_POLAPLIB_DIR}/run-polap-sh-create-depth-file.sh" "${_num1}" "${_num2}" "${_polap_var_mtcontigs_2_custom_depth_range}"
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
function _polap_seeds_create-automatic-depth-range {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	local _type=$1
	local -n result_ref=$2

	if [[ "${_arg_plastid}" == "on" ]]; then
		_polap_log2 "  using plastid depth range ..."
		local _command1="Rscript ${_POLAPLIB_DIR}/run-polap-r-plastid-determine-depth-range_${_type}.R \
				-t ${_polap_var_ga_annotation_all} \
				-c ${_polap_var_ga_annotation_cdf_table} \
				-o ${_polap_var_mtcontigs_1_custom_depth_range} \
        --plastid \
				2>$_polap_output_dest"
	else
		local _command1="Rscript ${_POLAPLIB_DIR}/run-polap-r-determine-depth-range_${_type}.R \
				-t ${_polap_var_ga_annotation_all} \
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
function _polap_seeds_preselect-contigs {

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
	_polap_log2 "    input2: ${_polap_var_ga_annotation_all}"
	_polap_log2 "    output: ${_polap_var_mtcontigs_preselection}"
	local _command1="Rscript ${_POLAPLIB_DIR}/run-polap-r-preselect-annotation.R  \
		--out ${_polap_var_mtcontigs_preselection} \
		--table ${_polap_var_ga_annotation_all} \
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
# step 3-1: creating GFA without sequence data
# step 3-2: extracting sequence part of GFA
# step 3-3: filtering GFA sequence part using depth range
################################################################################
function _polap_seeds_depthfilter-gfa {
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
		-S ${_polap_var_ga_contigger_edges_gfa} \
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
		_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/run-polap-r-depthfilter-gfa.R \
			--gfa ${_polap_var_mtcontigs_gfa_seq_part} \
			--depth ${_polap_var_mtcontigs_depth_range_graphfilter} \
			--out ${_polap_var_mtcontigs_gfa_seq_filtered} \
			2>$_polap_output_dest"
	else
		_polap_log0 "    no such file: ${_polap_var_mtcontigs_depth_range_graphfilter}"
		_polap_log0 "    ${_polap_var_mtcontigname} -> empty or no seed contigs"
		>"${_polap_var_mtcontigname}"
		_polap_log2_cat "${_polap_var_mtcontigname}"
		_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
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
function _polap_seeds_prepare-cc {
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
		${_polap_var_ga_contigger_edges_gfa} \
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
			sort | uniq >${_polap_var_mtcontigname}"
		_polap_log0 "  output1: ${_polap_var_mtcontigname}"
		_polap_log2_cat "${_polap_var_mtcontigname}"

		_polap_seeds_final-mtcontig

		_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
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
	_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/run-polap-r-prepare-cc.R \
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
function _polap_seeds_final-mtcontig {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# local _is_auto=$1 not using it anymore? for --length 1000000 option below

	_polap_log1 "  step 8-1: filtering by max length 1 Mb of contig seeds"
	_polap_log1 "    1 Mb length filtering for auto"
	_polap_log2 "    input1: ${_polap_var_ga_annotation_all}"
	_polap_log2 "    input2: ${_polap_var_mtcontigs_7mtcontigname}"
	_polap_log2 "    output1: ${_polap_var_mtcontig_table}"
	local _command1="Rscript ${_POLAPLIB_DIR}/run-polap-r-final-filter-mtcontig.R \
		-t ${_polap_var_ga_annotation_all} \
		-m ${_polap_var_mtcontigs_7mtcontigname} \
    -o ${_polap_var_mtcontig_table} \
    --length 1000000 \
		2>$_polap_output_dest"
	_polap_log3_pipe "${_command1}"

	_polap_seeds_final-seeds-mtcontig \
		"${_polap_var_mtcontigs_7mtcontigname}" \
		"${_polap_var_mtcontigs_annotation_table_seed}"
	# "${_polap_var_ga_annotation_depth_table_seed}"
	# _polap_log1 "  step 8-2: mt.contig.table with contig seed marks"
	# _polap_log2 "    input1: ${_polap_var_ga_annotation_all}"
	# _polap_log2 "    input2: ${_polap_var_ga_annotation_depth_table}"
	# _polap_log2 "    input3: ${_polap_var_mtcontigs_7mtcontigname}"
	# _polap_log2 "    output1: ${_polap_var_mtcontigs_annotation_table_seed}"
	# local _command1="Rscript ${_POLAPLIB_DIR}/run-polap-r-final-seed-mtcontig.R \
	# 	-t ${_polap_var_ga_annotation_all} \
	# 	-a ${_polap_var_ga_annotation_depth_table} \
	# 	-m ${_polap_var_mtcontigs_7mtcontigname} \
	#    -o ${_polap_var_mtcontigs_annotation_table_seed} \
	# 	2>$_polap_output_dest"
	# _polap_log3_pipe "${_command1}"

}

function _polap_seeds_final-seeds-mtcontig {
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"
	local _mt_contig_name=$1
	local _annotation_table_seed_target=$2

	_polap_log1 "  annotation table with contig seed marks"
	_polap_log2 "    input1: ${_polap_var_ga_annotation_all}"
	_polap_log2 "    input2: ${_polap_var_ga_annotation_depth_table}"
	_polap_log2 "    input3: ${_mt_contig_name}"
	_polap_log2 "    output1: ${_annotation_table_seed_target}"
	local _command1="Rscript ${_POLAPLIB_DIR}/run-polap-r-final-seed-mtcontig.R \
		-t ${_polap_var_ga_annotation_all} \
		-a ${_polap_var_ga_annotation_depth_table} \
		-m ${_mt_contig_name} \
    -o ${_annotation_table_seed_target} \
		2>$_polap_output_dest"
	_polap_log3_pipe "${_command1}"

}

function _polap_seeds_report-mtcontig {
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

################################################################################
# Prepare seed contigs using a Flye genome assembly.
################################################################################
function _run_polap_seeds { # select seed contigs
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
# Select contigs using seeds-graph
#
# Arguments:
#   -i ${_arg_inum}: source Flye (usually whole-genome) assembly number
#   -j ${_arg_jnum}: destination Flye organelle assembly number
#   --plastid
# Inputs:
#   ${_polap_var_ga_contigger_edges_gfa}
#   ${_polap_var_ga_annotation_all}
# Outputs:
#   ${_polap_var_mtcontigname} or more such files
# Menu:
#   bandage
#   annotation
# View:
#   <number> for the mt.contig.name-<number>
Example: $(basename "$0") ${_arg_menu[0]} -i 0 -j 1
Example: $(basename "$0") ${_arg_menu[0]} -i 0 -j 1 [--max-seeds ${_arg_max_seeds}]
Example: $(basename "$0") ${_arg_menu[0]} view -i 0 -j 2
Example: $(basename "$0") ${_arg_menu[0]} annotation
Example: $(basename "$0") ${_arg_menu[0]} bandage
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		if [[ "${_arg_markdown}" == "off" ]]; then
			_polap_log0 "$(basename $0) seeds annotation -i ${_arg_inum} -j ${_arg_jnum}"
			_polap_log0_column "${_polap_var_ga_annotation_depth_table_seed_target}"
		else
			csvtk space2tab "${_polap_var_ga_annotation_depth_table_seed_target}" | csvtk tab2csv | csvtk csv2md -a l,r,r,r,r,r,r >&3
		fi

		if [[ -s "${_polap_var_mtcontigname}" ]]; then
			_polap_log0 "---"
			if [[ "${_arg_log_stderr}" = "off" ]]; then
				paste -sd',' "${_polap_var_mtcontigname}" >&3
			else
				paste -sd',' "${_polap_var_mtcontigname}" >&2
			fi
		fi
		return $RETURN_SUCCESS
	fi

	# Backup mt.contig.name files.
	if [[ "${_arg_menu[1]}" == "backup" ]]; then
		# Format the current date as YYYY-MM-DD
		local folder_name=$(date +"%Y-%m-%d")
		# Create the folder with the formatted date as its name
		mkdir "${_polap_var_ga}/backup-mt.contig.name-$folder_name"
		mv "${_polap_var_ga}"/mt.contig.name-* "${_polap_var_ga}/backup-mt.contig.name-$folder_name"
		return $RETURN_SUCCESS
	fi

	if [[ "${_arg_menu[1]}" = "bandage" ]]; then
		# Prompt the user for input
		read -p "Enter edges using Bandage (e.g., edge_265, edge_520, edge_425): " edges

		# Convert the input string into new lines (replace ", " with newlines)
		local formatted_edges=$(echo "$edges" | tr ', ' '\n' | sed '/^ *$/d')
		echo "${formatted_edges}" >"${_polap_var_mtcontigname}"

		local file=$(basename ${_polap_var_mtcontigname})
		local number="${file#mt.contig.name-}"
		local _ga_annotation_depth_table_seed_target="${_polap_var_ga}/contig-annotation-depth-table-seed-${number}.txt"
		_polap_seeds_final-seeds-mtcontig \
			"${_polap_var_mtcontigname}" \
			"${_ga_annotation_depth_table_seed_target}"

		_polap_log0 "NEXT: $(basename "$0") assemble2 -o ${_arg_outdir} -i ${_arg_inum} -j ${_arg_jnum}"

		return $RETURN_SUCCESS
	fi

	if [[ "${_arg_menu[1]}" = "annotation" ]]; then
		# Prompt the user for input
		_polap_log0_cat "${_polap_var_mtcontigname}"
		_polap_seeds_final-seeds-mtcontig \
			"${_polap_var_mtcontigname}" \
			"${_polap_var_ga_annotation_depth_table_seed_target}"
		return $RETURN_SUCCESS
	fi

	# We initiate the process of selecting seed contigs.
	_polap_log0 "selecting seed contigs using the assembly graph: ${_arg_inum} (source) -> ${_arg_jnum} (target) ..."
	_polap_log1 "  input1: ${_polap_var_ga_contigger_edges_gfa}"
	_polap_log1 "  input2: ${_polap_var_ga_annotation_all}"

	rm -f "${_polap_var_ga}"/mt.contig.name-*

	# Create the array based on the value of _arg_plastid
	if [[ "$_arg_plastid" == "on" ]]; then
		local knum_array=(1 6)
	else
		local knum_array=({1..6}) # This creates an array of 1 through 6
	fi

	# 1. file_hashed.tmp
	# 2. unique_files_by_content.txt
	# 3. filtered_files.tmp

	_polap_log3_cmd mkdir -p "${_polap_var_ga_mtcontigs}"

	# Create a temporary file to store hashes and file paths
	local hash_file="${_polap_var_ga_mtcontigs}/file_hashes.tmp"
	>"$hash_file" # Ensure the file is empty
	for i in "${knum_array[@]}"; do
		_arg_knum=$i
		source "${_POLAPLIB_DIR}/polap-variables-mtcontigs.sh"
		_run_polap_seeds-graph

		if [[ -s "${_polap_var_mtcontigs_8mtcontigname}" ]]; then
			md5sum "${_polap_var_mtcontigs_8mtcontigname}" >>"$hash_file"
		fi
	done

	if [[ -s "${hash_file}" ]]; then
		_polap_log3_cat "${hash_file}"
	else
		return "${_POLAP_ERR_NO_SEEDS}"
	fi

	# This will give a list of unique files by content
	awk '!seen[$1]++' "$hash_file" >"${_polap_var_ga_mtcontigs}/unique_files_by_content.txt"

	# Temporary file to store files with fewer than 10 lines
	filtered_file="${_polap_var_ga_mtcontigs}/filtered_files.tmp"
	>"$filtered_file" # Ensure the file is empty
	while IFS=" " read -r hash file; do
		# Check if file is not empty and line count is less than _arg_max_seeds
		if [[ -f "$file" ]]; then
			local line_count=$(wc -l <"$file")
			if ((line_count < _arg_max_seeds)); then
				echo "$line_count $file" >>"$filtered_file"
			fi
		fi
	done <"${_polap_var_ga_mtcontigs}/unique_files_by_content.txt"

	# Initialize index counter
	local index=${_arg_jnum}

	# Loop through each unique file
	while IFS=" " read -r line_count file; do
		# Compute the new filename based on _arg_jnum and index
		local new_filename="${_polap_var_ga}/mt.contig.name-${index}"
		local new_filename_table="${_polap_var_ga}/contig-annotation-table-seed-${index}.txt"

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

	_polap_log0 "seed contig name files:"
	ls "${_polap_var_ga}"/mt.contig.name-* >&3

	for _mtcontigname in "${_polap_var_ga}"/mt.contig.name-*; do
		# Extract the <number> part using parameter expansion

		local file=$(basename ${_mtcontigname})
		local number="${file#mt.contig.name-}"
		local _ga_annotation_depth_table_seed_target="${_polap_var_ga}/contig-annotation-depth-table-seed-${number}.txt"

		_polap_seeds_final-seeds-mtcontig \
			"${_mtcontigname}" \
			"${_ga_annotation_depth_table_seed_target}"

		_polap_log0 "NEXT: $(basename "$0") assemble2 -o ${_arg_outdir} -i ${_arg_inum} -j ${number}"
	done

	_polap_log2 "Function end (${_arg_select_contig}): $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
