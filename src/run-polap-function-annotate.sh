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

source "$script_dir/run-polap-function-blast-genome.sh"
source "$script_dir/run-polap-function-count-gene.sh"

################################################################################
#
################################################################################
function _run_polap_edges-stats() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-ga.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Creates edges_stats.txt like contigs_stats.txt.
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
# Inputs:
#   ${_polap_var_contigger_gfa}
# Outputs:
#   ${_polap_var_contigger_edges_stats}
Example: $(basename $0) ${_arg_menu[0]} [-o ${ODIR}]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return
	fi

	mkdir -p "${_polap_var_ga_mtcontigs}"

	_polap_log1 "creating edges_stats.txt from graph_final.gfa ..."
	_polap_log2 "  creating GFA without sequence data: ${_polap_var_ga_gfa_all}"
	_polap_log2 "    input: ${_polap_var_contigger_gfa}"
	_polap_log2 "    output: ${_polap_var_ga_gfa_all}"
	gfatools view \
		-S "${_polap_var_contigger_gfa}" \
		>"${_polap_var_ga_gfa_all}" \
		2>"$_polap_output_dest"

	_polap_log2 "  extracting sequence part of GFA: ${_polap_var_ga_gfa_seq_part}"
	_polap_log2 "    input: ${_polap_var_ga_gfa_all}"
	_polap_log2 "    output: ${_polap_var_ga_gfa_seq_part}"
	grep "^S" "${_polap_var_ga_gfa_all}" >"${_polap_var_ga_gfa_seq_part}"

	# Filter edges in GFA using depths.
	_polap_log2 "  filtering GFA sequence part using depth range"
	_polap_log2 "    input1: ${_polap_var_ga_gfa_seq_part}"
	_polap_log2 "    output1: ${_polap_var_contigger_edges_stats}"
	Rscript "$script_dir"/run-polap-r-edges-stats.R \
		--gfa "${_polap_var_ga_gfa_seq_part}" \
		--out "${_polap_var_contigger_edges_stats}" \
		2>"$_polap_output_dest"

	if [[ -s "${_polap_var_contigger_edges_stats}" ]]; then
		_polap_log2_cat "${_polap_var_contigger_edges_stats}"
	fi

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
#
################################################################################
function _run_polap_depth-distribution() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-mtcontig.sh"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Draw distributions of annotation tables.
#
# Distribution by the cutoff of MT gene groups Copy average - 3 x SD
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
# Inputs:
#   ${_polap_var_annotation_table}
# Outputs:
#   ${_polap_var_annotation_table}.pdf
Example: $(basename $0) ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		_polap_log0_file "${_polap_var_annotation_table}".pdf
		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	Rscript "$script_dir"/run-polap-depth-distribution.R \
		-t "${_polap_var_annotation_table}" \
		-o "${_polap_var_annotation_table}".pdf \
		2>"$_polap_output_dest"

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Annotates the genome assembly.
# Defaults:
# Outputs:
################################################################################
function _run_polap_annotate() {
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
# Annotates the genome assembly.
#
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   $INUM
# Outputs:
Example: $(basename $0) ${_arg_menu[0]} [-i|--inum <arg>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		# local _polap_var_ga_annotation="${_polap_var_ga}/assembly_info_organelle_annotation_count.txt"
		# local _polap_var_ga_annotation_all="${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt"
		_polap_log0 column -t "${_polap_var_ga_annotation}"
		if [ -s "${_polap_var_ga_annotation}" ]; then
			column -t "${_polap_var_ga_annotation}" >&2
		fi

		_polap_log0 column -t "${_polap_var_ga_annotation_depth_table}"
		if [ -s "${_polap_var_ga_annotation_depth_table}" ]; then
			column -t "${_polap_var_ga_annotation_depth_table}" >&2
		fi

		_polap_log3_cmd touch "${ODIR}"
		_polap_log3_cmd ln -f -s "${_polap_var_contigger_gfa#*/}" "${ODIR}"/wga.gfa

		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "annotating contigs with mitochondrial and plastid genes on the assembly number ${INUM} ..."

	if [ -s "${_polap_var_ga_annotation_all}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log0 "  found: ${_polap_var_ga_annotation_all}, so skipping the annotation ..."
	else
		_run_polap_blast-genome
		_run_polap_count-gene
	fi

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Annotates edge sequences of a flye genome assembly.
# Defaults:
# Outputs:
################################################################################
function _run_polap_annotate-edge() {
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
# Annotates the genome assembly.
#
# Arguments:
#   -i $INUM: a Flye genome assembly number
# Inputs:
#   $INUM
# Outputs:
Example: $(basename $0) ${_arg_menu[0]} [-i|--inum <arg>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log0 column -t "${_polap_var_ga_annotation}"
		if [ -s "${_polap_var_ga_annotation}" ]; then
			column -t "${_polap_var_ga_annotation}" >&2
		fi

		_polap_log0 column -t "${_polap_var_ga_annotation_depth_table}"
		if [ -s "${_polap_var_ga_annotation_depth_table}" ]; then
			column -t "${_polap_var_ga_annotation_depth_table}" >&2
		fi

		_polap_log3_cmd touch "${ODIR}"
		_polap_log3_cmd ln -f -s "${_polap_var_contigger_gfa#*/}" "${ODIR}"/wga.gfa

		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "annotating edge sequence contigs with mitochondrial and plastid genes on the assembly number ${INUM} ..."

	if [ -s "${_polap_var_ga_annotation_all}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log0 "  found: ${_polap_var_ga_annotation_all}, so skipping the annotation ..."
	elif [ -s "${_polap_var_ga_annotation_all}" ] && [ "${_arg_redo}" = "on" ]; then
		_polap_log0 "  found: ${_polap_var_ga_annotation_all}, with --redo option, so backup the annotation ..."
		_polap_log3_cmd cp -p ${_polap_var_ga_annotation_all} ${_polap_var_ga_annotation_all_backup}
		_run_polap_blast-genome-edge
		_run_polap_count-gene-edge
	else
		_run_polap_blast-genome-edge
		_run_polap_count-gene-edge
	fi

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

function _run_polap_a() {
	_run_polap_annotate-edge
}

function _run_polap_dd() {
	_run_polap_depth-distribution
}