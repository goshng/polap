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
# Makes menu commands as empty files.
################################################################################
function _run_polap_make-menus() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	help_message=$(
		cat <<HEREDOC
# Makes menu files for easier command completion.
#
# Arguments: none
# Inputs: none
# Outputs: empty files with filenames of menus.
Example: $(basename $0) ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	_polap_log0 "creating empty menu files ..."

	cat "$script_dir"/*.sh |
		grep "^function _run_polap" |
		grep run_polap |
		grep -v run_polap_x |
		sed 's/function _run_polap_//' | sed 's/() {//' |
		parallel touch {}

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Makes menu commands as empty files.
# Creates menus prefixed with x.
################################################################################
function _run_polap_make-menus-all() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	help_message=$(
		cat <<HEREDOC
# Makes menu files including development version for easier command completion.
#
# Arguments: none
# Inputs: none
# Outputs: empty files with filenames of menus.
Example: $(basename $0) ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	_polap_log0 "creating empty menu files including development versions ..."

	cat "$script_dir"/*.sh |
		grep "^function _run_polap" |
		grep run_polap |
		sed 's/function _run_polap_//' | sed 's/() {//' |
		parallel touch {}

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Deletes menu commands.
# Leaves make-menus command.
################################################################################
function _run_polap_clean-menus() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	help_message=$(
		cat <<HEREDOC
# Clean up the menu files.
#
# Arguments: none
# Inputs: none
# Outputs: empty files with filenames of menus.
Example: $(basename $0) ${_arg_menu[0]}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	_polap_log0 "cleaning up the empty menu files ..."

	cat "$script_dir"/*.sh |
		grep "^function _run_polap" |
		grep run_polap |
		sed 's/function _run_polap_//' | sed 's/() {//' |
		parallel rm -f {}

	# Leave the make-menus empty file.
	touch make-menus

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Lists menu of POLAP.
# You need to execute make-menus menu if nothing is displayed.
################################################################################
function _run_polap_list() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	help_message=$(
		cat <<HEREDOC
# Lists menu of POLAP.
#
# You need to execute make-menus menu if nothing is displayed.
Example: $(basename "$0") ${_arg_menu[0]} [all|main|assemble1|annotate|assemble2|polish]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	if [[ ${_arg_menu[1]} == "all" ]]; then
		find . -maxdepth 1 -type f -empty -exec basename {} \; |
			sort >&2
	elif [[ ${_arg_menu[1]} == "main" ]]; then
		_polap_log0 assemble1
		_polap_log0 annotate
		_polap_log0 assemble2
		_polap_log0 flye-polishing
		_polap_log0 prepare-polishing
		_polap_log0 polish
	elif [[ ${_arg_menu[1]} == "assemble1" ]]; then
		_polap_log0 reset
		_polap_log0 summary-reads
		_polap_log0 total-length-long
		_polap_log0 find-genome-size
		_polap_log0 reduce-data
		_polap_log0 flye1
	elif [[ ${_arg_menu[1]} == "annotate" ]]; then
		_polap_log0 blast-genome
		_polap_log0 count-gene
	elif [[ ${_arg_menu[1]} == "assemble2" ]]; then
		_polap_log0 select-reads
		_polap_log0 flye2
	elif [[ ${_arg_menu[1]} == "polish" ]]; then
		_polap_log0 flye-polishing
		_polap_log0 prepare-polishing
		_polap_log0 polish
	else
		_polap_log0 "${help_message}"
	fi

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
