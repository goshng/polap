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
# Bandage subcommand creates png figure files for the given assembly graphs.
# It uses the Bandage tool to generate images from assembly graphs.
#
# Functions:
#   _run_polap_bandage
#
# Usage:
#   polap bandage png input.gfa output.png [label|random]
#
# TEST-SCC: not yet.
# TEST-DOC: not yet.
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

function _run_polap_bandage {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	# Print help message if requested
	help_message=$(
		cat <<HEREDOC
# Use Bandage command to generate a genome assembly graph
#
# Inputs:
#  input.gfa: Input assembly graph file in GFA format.
#
# Outputs:
#  output.png: PNG image file generated from the assembly graph.
#
# Menu:
#   png: Generate a PNG image from the assembly graph.
#
# See:
Example: $(basename $0) ${_arg_menu[0]} png input.gfa output.png [label|random]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "png" ]]; then
		local _infile="${_arg_menu[2]}"
		local _outfile="${_arg_menu[3]}"
		local _color="${_arg_menu[4]}"
		if [[ "${_color}" == "fourthfile" ]]; then
			# gray color
			Bandage image "${_infile}" \
				"${_outfile}" \
				--colour uniform \
				--unicolpos \#EEEEEE \
				--singlearr \
				--toutline 1
		elif [[ "${_color}" == "label" ]]; then
			# gray color
			Bandage image "${_infile}" \
				"${_outfile}" \
				--colour uniform \
				--unicolpos \#EEEEEE \
				--singlearr \
				--toutline 1 \
				--names \
				--lengths \
				--depth \
				--fontsize 3
		else
			# random color
			Bandage image "${_infile}" \
				"${_outfile}" \
				--colour random \
				--singlearr \
				--toutline 1 \
				--names \
				--lengths \
				--depth \
				--fontsize 3
		fi
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
