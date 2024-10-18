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
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

################################################################################
# Counts genes annotated on a genome assembly.
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   ${_polap_var_ga}/30-contigger/contigs_stats.txt
#   $MTGENECOUNT
#   $PTGENECOUNT
# Outputs:
#   ${_polap_var_ga}/mt.contig.name-1
#   ${_polap_var_ga}/mt.contig.name-2
#   ${_polap_var_ga}/assembly_info_organelle_annotation_count.txt
#   ${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt
#   ${_polap_var_ga}/contig-annotation-table.txt
################################################################################
function _run_polap_count-gene() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-ga.sh"

	help_message=$(
		cat <<HEREDOC
# Counts genes annotated on a genome assembly.
#
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   ${_polap_var_ga}/30-contigger/contigs_stats.txt
#   ${_polap_var_ann_MTGENECOUNT}
#   ${_polap_var_ann_PTGENECOUNT}
# Outputs:
#   ${_polap_var_ga}/mt.contig.name-1
#   ${_polap_var_ga}/mt.contig.name-2
#   ${_polap_var_ga}/assembly_info_organelle_annotation_count.txt
#   ${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt
#   ${_polap_var_ga}/contig-annotation-table.txt 
#   ${_polap_var_ga}/contig-annotation-depth-table.txt 
# Usage:
#   assembly graph: ${_polap_var_ga}"/30-contigger/graph_final.gfa
# 	column -t ${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt | less -S
#   column -t ${_polap_var_ga}/assembly_info_organelle_annotation_count.txt | less -S
#   column -t ${_polap_var_ga}/contig-annotation-table.txt | less -S
# Example file: ${_polap_var_ga}/mt.contig.name-1
# 	edge_1
#   edge_2
#   edge_3
# edit ${_polap_var_ga}/mt.contig.name-1 for mtDNA contig candidates
# edit ${_polap_var_ga}/mt.contig.name-<destination flye number> for mtDNA contig candidates
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		# Display the BLAST genome output.

		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "counting mitochondrial and plastid genes on the assembly number ${INUM} ..."

	_polap_log0 "  input1: ${_polap_var_ga}/30-contigger/contigs_stats.txt"
	_polap_log0 "  input2: ${_polap_var_ann_MTGENECOUNT}"
	_polap_log0 "  input3: ${_polap_var_ann_PTGENECOUNT}"

	Rscript "$script_dir"/run-polap-mtcontig.R \
		"${_polap_var_ga}" \
		"${_polap_var_ga}"/mt.contig.name \
		"${_polap_var_ga_annotation}" \
		"${_polap_var_ga_annotation_all}" \
		--contigger \
		>${_polap_output_dest} 2>&1

	_polap_log0 "  output1: ${_polap_var_ga_annotation}"
	_polap_log0 "  output2: ${_polap_var_ga_annotation_all}"
	_polap_log0 "  output3: ${_polap_var_ga_annotation_table}"
	_polap_log0 "  output4: ${_polap_var_ga_annotation_depth_table}"

	if [[ ${_arg_test} == "on" ]]; then
		_polap_log0 "creating ${_polap_var_ga}/mt.contig.name-1 for testing purpose ..."
		_polap_log0 "  you would have to edit it for a real data-set."
		echo edge_1 >"${_polap_var_ga}"/mt.contig.name-1
		echo edge_2 >>"${_polap_var_ga}"/mt.contig.name-1
		echo edge_3 >>"${_polap_var_ga}"/mt.contig.name-1
	fi

	ANUMNEXT=$((INUM + 1))
	echoerr NEXT: $(basename "$0") select-reads -o "$ODIR" [-i $INUM] [-j $ANUMNEXT]
	echoerr NEXT: $(basename "$0") assemble2 -o "$ODIR" [-i $INUM] [-j $ANUMNEXT]

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

function _run_polap_cg() {
	_run_polap_count-gene
}
################################################################################
# Counts genes annotated on a genome assembly.
#
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   ${_polap_var_ga}/30-contigger/edges_stats.txt
#   $MTGENECOUNT
#   $PTGENECOUNT
# Outputs:
#   ${_polap_var_ga}/mt.contig.name-1
#   ${_polap_var_ga}/mt.contig.name-2
#   ${_polap_var_ga}/assembly_info_organelle_annotation_count.txt
#   ${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt
#   ${_polap_var_ga}/contig-annotation-table.txt
################################################################################
function _run_polap_count-gene-edge() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-ga.sh"

	help_message=$(
		cat <<HEREDOC
# Counts genes annotated on a genome assembly.
#
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   ${_polap_var_ga}/30-contigger/edges_stats.txt
#   ${_polap_var_ann_MTGENECOUNT}
#   ${_polap_var_ann_PTGENECOUNT}
# Outputs:
#   ${_polap_var_ga}/mt.contig.name-1
#   ${_polap_var_ga}/mt.contig.name-2
#   ${_polap_var_ga}/assembly_info_organelle_annotation_count.txt
#   ${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt
#   ${_polap_var_ga}/contig-annotation-table.txt 
#   ${_polap_var_ga}/contig-annotation-depth-table.txt 
# Usage:
#   assembly graph: ${_polap_var_ga}"/30-contigger/graph_final.gfa
# 	column -t ${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt | less -S
#   column -t ${_polap_var_ga}/assembly_info_organelle_annotation_count.txt | less -S
#   column -t ${_polap_var_ga}/contig-annotation-table.txt | less -S
# Example file: ${_polap_var_ga}/mt.contig.name-1
# 	edge_1
#   edge_2
#   edge_3
# edit ${_polap_var_ga}/mt.contig.name-1 for mtDNA contig candidates
# edit ${_polap_var_ga}/mt.contig.name-<destination flye number> for mtDNA contig candidates
Example: $(basename "$0") ${_arg_menu[0]} [-i|--inum <arg>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		# Display the BLAST genome output.

		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "counting mitochondrial and plastid genes on the assembly number ${INUM} ..."

	_polap_log0 "  input1: ${_polap_var_ga}/30-contigger/edges_stats.txt"
	_polap_log0 "  input2: ${_polap_var_ann_MTGENECOUNT}"
	_polap_log0 "  input3: ${_polap_var_ann_PTGENECOUNT}"
	_polap_log0 "  output1: ${_polap_var_ga_annotation}"
	_polap_log0 "  output2: ${_polap_var_ga_annotation_all}"
	_polap_log0 "  output3: ${_polap_var_ga_annotation_table}"
	_polap_log0 "  output4: ${_polap_var_ga_annotation_depth_table}"

	_polap_log3_pipe "Rscript $script_dir/run-polap-r-mtcontig.R \
		--flyeout-edges-stats ${_polap_var_contigger_edges_stats} \
    --mt-gene-count ${_polap_var_ann_MTGENECOUNT} \
    --pt-gene-count ${_polap_var_ann_PTGENECOUNT} \
		--out-annotation ${_polap_var_ga_annotation} \
		--out-annotation-all ${_polap_var_ga_annotation_all} \
		--out-annotation-table ${_polap_var_ga_annotation_table} \
		--out-annotation-depth-table ${_polap_var_ga_annotation_depth_table} \
		--contigger \
		2>${_polap_output_dest}"

	if [[ -s "${_polap_var_ga_annotation_all}" ]]; then
		_polap_log3_cat "${_polap_var_ga_annotation_all}"
	fi

	if [[ ${_arg_test} == "on" ]]; then
		_polap_log0 "creating ${_polap_var_ga}/mt.contig.name-1 for testing purpose ..."
		_polap_log0 "  you would have to edit it for a real data-set."
		echo edge_1 >"${_polap_var_ga}"/mt.contig.name-1
		echo edge_2 >>"${_polap_var_ga}"/mt.contig.name-1
		echo edge_3 >>"${_polap_var_ga}"/mt.contig.name-1
	fi

	ANUMNEXT=$((INUM + 1))
	echoerr NEXT: $(basename "$0") select-reads -o "$ODIR" [-i $INUM] [-j $ANUMNEXT]
	echoerr NEXT: $(basename "$0") assemble2 -o "$ODIR" [-i $INUM] [-j $ANUMNEXT]

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

function _run_polap_cg() {
	_run_polap_count-gene
}
