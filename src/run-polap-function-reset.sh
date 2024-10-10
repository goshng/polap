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
# Template for an external shell script
#
# You could use this function template to create a new menu.
# Rename template and delete x in the name. You could execute such menu
# but such menus are not created as empty files by make-menus menu.
################################################################################

################################################################################
# Ensure that the current script is sourced only once
source "$script_dir/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

source "$script_dir/run-polap-function-utilities.sh"
source "$script_dir/run-polap-function-menus.sh"

################################################################################
# Initializes polap analysis in a starting folder,
# creating an output folder.
# Arguments:
#   -o $ODIR
# Inputs: nothing
# Outputs:
#   $ODIR
################################################################################
function _run_polap_reset() {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	help_message=$(
		cat <<HEREDOC
# Initializes polap analysis in a starting folder by creating the output folder.
#
# Arguments:
#   -o ${ODIR}: the output folder
# Inputs: none
# Outputs:
#   ${ODIR}
Example: $(basename "$0") ${_arg_menu[0]} [-o|--outdir ${ODIR}]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	if ! run_check1; then
		error_polap_conda
		exit $EXIT_ERROR
	fi

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [[ -d "${ODIR}" ]]; then
			ls "${ODIR}" >&2
		else
			_polap_log0 "No such output folder: ${ODIR}"
		fi
		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	if [ -d "$ODIR" -a "${_arg_yes}" = "off" ]; then
		while true; do
			read -p "Folder [$ODIR] already exists. Do you want to delete it? [y/n] " yn
			case $yn in
			[Yy]*)
				rm -rf "$ODIR"
				break
				;;
			[Nn]*)
				exit $EXIT_FAIL
				;;
			*) _polap_log0 "Please answer yes or no." ;;
			esac
		done
	else
		rm -rf "$ODIR"
	fi

	mkdir -p "$ODIR"
	_polap_log0 "creating output folder [$ODIR] ..."
	_polap_log1 "  Your output folder [$ODIR] is created."
	if [ "$ODIR" != "o" ]; then
		_polap_log1 "  Use -o $ODIR option in all subsequent analysis"
		_polap_log1 "  because your output folder is not the default of 'o'."
	fi
	_run_polap_make-menus

	_polap_log1 NEXT: $(basename "$0") total-length-long -o "$ODIR" -l ${_arg_long_reads}
	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
