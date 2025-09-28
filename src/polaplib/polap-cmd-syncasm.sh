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

function _run_polap_syncasm {
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
  Get organelle genome sequences:
    polap ${polap_cmd} -l l.fq

TODO:
  Dev.

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	if [[ "${_arg_steps_is}" == "off" ]]; then
		_arg_steps_include="all"
	fi

	_polap_log0 "oatk-c: ${_arg_oatk_c}"
	_polap_log0 "parallel: ${_arg_parallel}"

	_polap_lib_conda-ensure_conda_env polap-oatk || exit 1
	_polap_lib_syncasm -v \
		-l "${_arg_long_reads}" \
		-o "${_arg_outdir}/syncasm1" \
		--steps "${_arg_steps_include}" \
		-t "${_arg_threads}" \
		--k "${_arg_oatk_k}" \
		--s "${_arg_oatk_s}" \
		--c "${_arg_oatk_c}" \
		--a "${_arg_oatk_a}" \
		--max-bubble "${_arg_oatk_max_bubble}" \
		--max-tip "${_arg_oatk_max_tip}" \
		--weak-cross "${_arg_oatk_weak_cross}" \
		--unzip-round "${_arg_oatk_unzip_round}" \
		--oatkdb "${_arg_oatk_oatkdb}" \
		--max-total-bp "3000000" \
		--min-seq-len "10000" \
		--minimap2-min-identity 0.75 \
		--minimap2-min-qcov 0.50 \
		--read-min-qcov 0.80 \
		$([[ "$_arg_plastid" == "true" ]] && echo "--plastid" || echo "") \
		$([[ "$_arg_oatk_no_read_ec" == "true" ]] && echo "--no-read-ec" || echo "") \
		$([[ "$_arg_oatk_no_hpc" == "true" ]] && echo "--no-hpc" || echo "") \
		$([[ "$_arg_parallel" == "true" ]] && echo "--parallel" || echo "")
	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
