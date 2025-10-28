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

function _run_polap_template {
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
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	_polap_lib_help-maybe-show3 "$polap_cmd" help_message || return 0
	# if [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]]; then
	# 	local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]}")
	# 	man "$manfile" >&3
	# 	rm -f "$manfile"
	# 	return
	# fi

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	_polap_lib_conda-ensure_conda_env polap || exit 1

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

	local OUTDIR="${_arg_outdir}/polish-longshort"
	_polap_log3_cmd rm -rf "$OUTDIR"

	_polap_log3_cmd bash "${_POLAPLIB_DIR}/polap-bash-polish-hybrid-racon-fmlrc2-polypolish.sh" \
		--ont "${_arg_long_reads}" \
		--sr1 "${_arg_short_read1}" \
		--sr2 "${_arg_short_read2}" \
		--fasta "${_arg_infile}" \
		--outdir "${OUTDIR}" \
		--threads "${_arg_half_threads}" \
		--rounds 2 --min-ident 0.40 --min-alen 2000 \
		${_arg_verbose_str}

	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
