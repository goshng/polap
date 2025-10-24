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

_polap_var_function_verbose=4

# Each _run_polap_ function has these 3 lines to set the verbose level.
# This source is intended to be used as a global so that each _run_polap_
# function does not need to have these 3 lines.
# However, it is not used much.
#
# Set verbosity level: stderr if verbose >= 2, otherwise discard output
_polap_output_dest="/dev/null"
# echo "verbose1: $_arg_verbose"
if [ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ]; then
	_polap_output_dest="/dev/stderr"
fi
# echo "verbose2: $_arg_verbose"

# Version: v0.2.2
trace_functions_fline() {
	local outfd="${1:-2}"

	# Ensure traps run inside functions/command substitutions/subshells
	# Remember prior state to restore later
	local prev_functrace
	prev_functrace="$(set +o | grep -E '^functrace')"
	set -o functrace # or: set -T

	: "${TRACE_FN_DEPTH:=0}"

	trap '{
    local depth=${#FUNCNAME[@]}
    local prev=${TRACE_FN_DEPTH:-0}

    if (( depth > prev )); then
      local fn="${FUNCNAME[0]:-MAIN}"
      [[ $fn == trace_functions_fline || $fn == untrace_functions ]] || {
        local call_file="${BASH_SOURCE[1]:-${BASH_SOURCE[0]:-?}}"
        local call_line="${BASH_LINENO[0]:-${LINENO}}"
        local indent; printf -v indent "%*s" $((depth-1)) ""
        builtin printf "%s→ %s (%s:%s)\n" "${indent// /·}" "$fn" "$call_file" "$call_line" >&'"$outfd"'
      }
    fi
    TRACE_FN_DEPTH=$depth
  }' DEBUG

	trap '{
    local fn="${FUNCNAME[0]:-MAIN}"
    [[ $fn == trace_functions_fline || $fn == untrace_functions ]] && return
    local depth=${#FUNCNAME[@]}
    local indent; printf -v indent "%*s" $((depth-1)) ""
    builtin printf "%s← %s\n" "${indent// /·}" "$fn" >&'"$outfd"'
  }' RETURN

	# Store how to restore functrace
	TRACE_FN_PREV_FUNCTRACE="$prev_functrace"
}

untrace_functions() {
	trap - DEBUG RETURN
	# Restore previous functrace setting
	if [[ -n "${TRACE_FN_PREV_FUNCTRACE:-}" ]]; then
		eval "$TRACE_FN_PREV_FUNCTRACE" # e.g., "set +o functrace"
		unset TRACE_FN_PREV_FUNCTRACE
	fi
	unset TRACE_FN_DEPTH
}
