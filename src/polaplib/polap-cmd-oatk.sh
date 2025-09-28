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

function _run_polap_oatk {
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
  polap ${polap_cmd} - assemble organelle genomes using Oatk's components.

Synopsis:
  polap ${polap_cmd} [options]

Description:
  while Oatk (and SyncAsm) was built around HiFi, you can use its components on ONT successfully if you adapt a few things. The “HiFi-only” idea comes from early assumptions about error rates and graph heuristics; once you account for ONT’s profile, the same closed-syncmer/sparse-DBG machinery is perfectly usable.

Options:
  -l FASTQ
    reads data file

Examples:
  Get organelle genome sequences:
    polap ${polap_cmd} -l l.fq

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

	# read the long read data
	# annotate reads for mt or pt anchor
	# call polap-bash-autotune-oatk.sh: --wgs-mode or --no-wgs-mode
	# We want to use oatk's components to assemble organelle genomes.

	# base directory
	local oatkdir="${_arg_outdir}/oatk"
	rm -rf "${oatkdir}"
	mkdir -p "${oatkdir}"

	# take the input long-read data and save it in outdir for use
	# _polap_lib_fastq-normalize-filename "${_arg_long_reads}"
	# _arg_long_reads="${_arg_outdir}/tmp/l.fq"
	# if [[ ! -s "${_arg_long_reads}" ]]; then
	# 	_polap_log0 "[ASSERT] no file: ${_arg_long_reads}"
	# 	return
	# fi

	# annotate reads for mt/pt anchors
	_polap_lib_readassemble-annotate-read-pt "oatk"
	_polap_lib_readassemble-get-annotated-read "oatk"

	################################################################################
	# local CUT_AUTO="${_POLAPLIB_DIR}/polap-py-syncfilter-cut-auto.py"
	# local PPR_SEL="${_POLAPLIB_DIR}/polap-py-syncmer-connectivity-select-mt.py"
	# local ENSEMBLE_PY="${_POLAPLIB_DIR}/polap-py-ensemble-mt.py"
	# local REPORT_PY="${_POLAPLIB_DIR}/polap-py-ensemble-report.py"
	#
	# polap-bash-select-mt.sh
	#
	# polap-bash-oatk-tune-c.sh
	# polap-bash-oatk-score-seeds.sh
	#
	# polap-bash-autotune-oatk.sh
	# <- polap-bash-oatk-tune-c.sh
	# <- polap-bash-oatk-score-seeds.sh
	# <- polap-bash-select-mt.sh
	#
	# polap-launch-autotune.sh <- polap-bash-autotune-oatk.sh
	#
	# How do I call polap-launch-autotune.sh?

	_polap_lib_conda-ensure_conda_env polap-oatk || exit 1

	local annotatedir="${_arg_outdir}/annotate-read-oatk"
	bash "${_POLAPLIB_DIR}/polap-launch-autotune.sh" \
		--preset hifi \
		--label mt \
		--verbose \
		--reads "${_arg_long_reads}" \
		--mt-anchors "${annotatedir}/mt.id.all.txt" \
		--pt-anchors "${annotatedir}/pt.id.all.txt" \
		--out-prefix "${_arg_outdir}/oatk/mt"

	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
