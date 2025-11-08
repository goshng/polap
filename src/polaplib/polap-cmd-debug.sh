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
# Functions for subcommand template ...
# Describe what they are and what they do.
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

function _run_polap_debug {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	local polap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  polap ${polap_cmd} - annotate rougly reads with organelle genes

Synopsis:
  polap ${polap_cmd} [options]

Description:
  polap ${polap_cmd} uses plastid and organelle genes to annotate reads
  using minimap2.

Options:
  -l FASTQ
    reads data file

Examples:
  Test fail with bash:
    polap ${polap_cmd} --lang bash --mode cmd

    polap ${polap_cmd} --lang bash --mode exit

    polap ${polap_cmd} --lang bash --mode nounset

    polap ${polap_cmd} --lang bash-py

    polap ${polap_cmd} --lang bash-r

    polap ${polap_cmd} --lang r

    polap ${polap_cmd} --lang py

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	_polap_lib_help-maybe-show3 "$polap_cmd" help_message || return 0

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	local OUTDIR="${_arg_outdir}/debug"
	_polap_log3_cmd rm -rf "$OUTDIR"

	_polap_lib_conda-ensure_conda_env polap || exit 1

	if [[ "${_arg_lang}" == "bash" ]]; then
		bash "${_POLAPLIB_DIR}/polap-bash-test-crash-bash.sh" --mode "${_arg_mode}" --code 7
	elif [[ "${_arg_lang}" == "bash-py" || "${_arg_lang}" == "bash-python" ]]; then
		bash "${_POLAPLIB_DIR}/polap-bash-test-crash-python.sh" --code 7
	elif [[ "${_arg_lang}" == "bash-r" || "${_arg_lang}" == "bash-R" ]]; then
		bash "${_POLAPLIB_DIR}/polap-bash-test-crash-r.sh" --code 7
	elif [[ "${_arg_lang}" == "py" || "${_arg_lang}" == "python" ]]; then
		python3 "${_POLAPLIB_DIR}/polap-py-test-crash.py" --code 7
	elif [[ "${_arg_lang}" == "r" || "${_arg_lang}" == "R" ]]; then
		Rscript --vanilla "${_POLAPLIB_DIR}/polap-r-test-crash.R" --code 7
	fi

	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
