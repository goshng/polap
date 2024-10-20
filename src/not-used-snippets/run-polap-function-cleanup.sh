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
# Clean up the ${ODIR}.
################################################################################
function _run_polap_cleanup() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	local FDIR="$ODIR"/$INUM
	local MTCONTIGNAME="$FDIR"/mt.contig.name-"$JNUM"
	local _polap_var_mtcontigs="$FDIR"/"$JNUM"/mtcontigs
	local _polap_var_assembly_graph_final_gfa="${FDIR}/30-contigger/graph_final.gfa"
	local _polap_var_annotation_table="${FDIR}/assembly_info_organelle_annotation_count-all.txt"

	local _polap_var_mtcontig_base="${_polap_var_mtcontigs}/1-mtcontig"
	local _polap_var_mtcontig_stats="${_polap_var_mtcontigs}/1-mtcontig.stats.txt"
	local _polap_var_mtcontig_annotated="${_polap_var_mtcontigs}/1-mtcontig.annotated.txt"

	local _polap_var_mtcontigs_mt_stats="${_polap_var_mtcontigs}/1-mtcontig.mt.stats.txt"
	local _polap_var_mtcontigs_pt_stats="${_polap_var_mtcontigs}/1-mtcontig.pt.stats.txt"

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Clean up the ${ODIR}.
#
# Arguments:
#   -o $ODIR: the output directory
#   -o $ODIR: the output directory
#
# Inputs:
#   ${ODIR}
#
# Outputs:
#
Example: $(basename $0) ${_arg_menu[0]} [-i|--inum <arg>] [-j|--jnum <arg>] [--select-contig <number>]
Example: $(basename $0) ${_arg_menu[0]} -o PRJNA914763 -i 0 -j 5 --select-contig 5
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
