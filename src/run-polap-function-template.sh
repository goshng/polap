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

function _run_polap_template() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "$script_dir/polap-variables-bioproject.sh" # '.' means 'source'

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Template for an external shell script
#
# Arguments:
#   -i $INUM: source Flye (usually whole-genome) assembly number
#
# Inputs:
#   ${_polap_var_annotation_table}
#
# Outputs:
#   ${MTCONTIGNAME}
#
# See:
#   run-polap-select-contigs-by-table-1.R for the description of --select-contig option
Example: $(basename $0) ${_arg_menu[0]} [-i|--inum <arg>] [-j|--jnum <arg>] [--select-contig <number>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	echo "verbose level: ${_arg_verbose}" >&2
	echoall "command: $0"
	echoall "function: $FUNCNAME"
	echoall "menu2: [$1]"
	echoall "menu3: [$2]"
	echoerr "LOG: echoerr"
	echoall "LOG: echoall"

	echoerr "LOG: echoerr"
	verbose_echo 0 "Log level   - screen        polap.log file" 1>&2
	_polap_log0 "Log level 0 - nothing        minimal log - --quiet"
	_polap_log1 "Log level 1 - minimal        step info and main io files"
	_polap_log2 "Log level 2 - main io files  inside of function: file input/output --verbose"
	_polap_log3 "Log level 3 - files inside   all log or details of file contents --verbose --verbose"
	_polap_log0_file "log0.file: main assembly input/output"
	_polap_log1_file "log1.file: step main input/output"
	_polap_log2_file "log2.file: inside detail input/output"
	_polap_log3_file "log3.file: all input/output"

	_polap_log0 "var: ${_polap_var_apple}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
