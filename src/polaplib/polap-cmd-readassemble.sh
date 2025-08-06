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

function _run_polap_readassemble {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	help_message=$(
		cat <<'EOF'
Name:
  polap readassemble - annotate reads before organelle genome assembly

Synopsis:
  polap readassemble [options]

Description:
  polap readassemble uses organelle genes to annotate reads before organelle genome assembly.

Options:
  -l FASTQ
    long reads data file

  --plastid
    assembly mtDNA instead of mtDNA

  --animal
    assemble animal mtDNA

  --nano-raw [default]
  --pacbio-hifi

Examples:
  Assemble plant mitochondrial sequences:
    polap readassemble -l l.fq

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

	# Three cases
	# 1. plastid
	# 2. mitochondrial or mitochondrial + noncoding
	# 3. animial mitochondrial
	if [[ "${_arg_plastid}" == "on" ]]; then
		_polap_log1 "Read-assemble ptDNA"
		_polap_readassemble-pt
	elif [[ "${_arg_animal}" == "on" ]]; then
		_polap_log1 "Read-assemble animal mtDNA"
		_polap_readassemble-mt
	else
		if [[ "${_arg_noncoding}" == "on" ]]; then
			_polap_log1 "Read-assemble plant mtDNA with mitochondrial noncoding regions"
			_polap_readassemble-nt
		else
			_polap_log1 "Read-assemble plant mtDNA without mitochondrial noncoding regions"
		fi
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

_polap_readassemble-nt() {
	_polap_lib_readassemble-annotate-read-nt
	_polap_lib_readassemble-assemble-annotated-read-nt
}

_polap_readassemble-pt() {
	_polap_lib_readassemble-annotate-read-pt
	_polap_lib_readassemble-assemble-annotated-read-pt
}

_polap_readassemble-mt() {
	_polap_lib_readassemble-annotate-read-mt
	_polap_lib_readassemble-assemble-annotated-read-mt
}
