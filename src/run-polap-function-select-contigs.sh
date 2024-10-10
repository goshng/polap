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

source "$script_dir/run-polap-function-select-contigs-by.sh"

################################################################################
# Select seed contigs using multiple methods
################################################################################
function _run_polap_select-contigs() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	help_message=$(
		cat <<HEREDOC
# Select seed contigs using multiple methods.
#
# Arguments:
#   -o $ODIR
#   -i $INUM: source Flye (usually whole-genome) assembly number
# Inputs:
#   $ODIR/$INUM
# Outputs:
#   $ODIR/$INUM/1
#   $ODIR/$INUM/2
#   $ODIR/$INUM/3
#   $ODIR/$INUM/4
#   $ODIR/$INUM/5
Example: $(basename $0) ${_arg_menu[0]} [-o $ODIR] [-i <number>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [[ "${_arg_menu[2]}" == "simple" ]]; then
			wc "${ODIR}/${INUM}"/mt.contig.name-* >&2
			exit $EXIT_SUCCESS
		fi

		for i in "${_arg_select_contig_numbers[@]}"; do

			local MTCONTIGNAME="${ODIR}/${INUM}/mt.contig.name-${i}"

			if [ -s "$MTCONTIGNAME" ]; then
				_arg_select_contig="${i}"
				JNUM="${i}"
				_run_polap_select-contigs-by
			fi
		done

		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "selecting contigs from the genome assembly number ${INUM} ..."

	# Loop over numbers from 1 to 5
	for i in "${_arg_select_contig_numbers[@]}"; do
		# Call the function corresponding to the current number (index is i-1)
		if [[ "${INUM}" -ne 0 ]]; then
			_polap_log0 "ASSERT: only -i 0 is implemented"
			exit $EXIT_FAIL
		fi
		local MTCONTIGNAME="${ODIR}/${INUM}/mt.contig.name-${i}"

		if [ -e "$MTCONTIGNAME" ] && [ "${_arg_redo}" = "off" ]; then
			_polap_log0 "  found: ${MTCONTIGNAME}, so skipping the select-contig step for ${MTCONTIGNAME}"
		else
			_arg_select_contig="${i}"
			JNUM="${i}"
			_run_polap_select-contigs-by
		fi

		# check the mt.contig.name-1
		if [ -s "$MTCONTIGNAME" ]; then
			_polap_log1_file "${MTCONTIGNAME}"
		else
			_polap_log0 "  $MTCONTIGNAME is empty; choose seed contigs by yourself."
		fi

	done

	if [ "${_arg_verbose}" -ge 0 ]; then
		wc "${ODIR}/${INUM}"/mt.contig.name-* >&2
	fi

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

function _run_polap_sc() {
	_run_polap_select-contigs
}
