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

function _run_polap_miniassemble {
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
  polap ${polap_cmd} - assemble organelle genomes using miniasm as a reference generator

Synopsis:
  polap ${polap_cmd} [options]

Description:
  A test replicate of readassemble.

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

	_polap_lib_conda-ensure_conda_env polap || exit 1

	_arg_steps_include="4,5"

	local _include="${_arg_steps_include}"
	local _exclude="${_arg_steps_exclude}" # Optional range or list of steps to exclude
	local _step_array=()

	_step_array=($(_polap_parse_steps "${_include}" "${_exclude}"))
	_rstatus="$?"
	if [[ "${_rstatus}" -ne 0 ]]; then
		_polap_log0 "ERROR: Error parsing steps."
		return "${_POLAP_ERR_CMD_OPTION_STEPS}"
	fi

	if _polap_contains_step 1 "${_step_array[@]}"; then
		_polap_log0 "step 1"
		# assemble pt
		_arg_plastid="on"
		if [[ "${_arg_readassemble_pt}" == "on" ]]; then
			_polap_readassemble-pt
		else
			_arg_long_reads="${_arg_outdir}/ld.fq"
		fi
	fi

	if _polap_contains_step 2 "${_step_array[@]}"; then
		_polap_log0 "step 2"
		# select mt reads for a very early seeds
		_arg_long_reads="${_arg_long_reads_original}"
		_arg_plastid="off"
		_polap_lib_readassemble-select-organelle-reads mtseed
		_polap_lib_file-cleanup -d "${_arg_outdir}/annotate-read-mtseed" -s 5M -a rm
	fi

	if _polap_contains_step 3 "${_step_array[@]}"; then
		_polap_log0 "step 3"
		# generate mt seed contigs
		_arg_pt_ref="${_arg_outdir}/pt.1.fa"
		_polap_assert '[[ -s "${_arg_pt_ref}" ]]' "pt ref must exist, '${_arg_pt_ref}'"
		rm -rf "${_arg_outdir}/mtseed/mt"{1..9}

		bash "${_POLAPLIB_DIR}/polap-bash-fast-mtseed-ont.sh" \
			-r "${_arg_long_reads}" \
			-o "${_arg_outdir}/mtseed" \
			-p "${_arg_pt_ref}" \
			--pt-origin "${_arg_outdir}/annotate-read-mtseed/pt.id.all.txt" \
			--mt-origin "${_arg_outdir}/annotate-read-mtseed/mt.id.all.txt" \
			-n "busco_downloads/lineages/viridiplantae_odb12/refseq_db.faa.gz" \
			-t "${_arg_threads}" \
			--use-parallel \
			--no-do-polap \
			${_arg_verbose_str} \
			--step "${_arg_steps_include}"
	fi

	# _polap_lib_file-cleanup -d "${_arg_outdir}/mtseed" -s 5M -a rm

	# assemble mtDNA using the miniasm seeds

	# assemble mtDNA using the flye seeds
	if _polap_contains_step 4 "${_step_array[@]}"; then
		_polap_log0 "step 4"
		local FDIR_NAME="07-flye"
		local FDIR="${_arg_outdir}/$FDIR_NAME"

		local _backup_outdir="${_arg_outdir}"
		_arg_inum="${FDIR_NAME}"
		_arg_outdir="${_arg_outdir}/mtseed"
		_polap_lib_readassemble-miniasm mtseed
		_arg_outdir="$_backup_outdir"
	fi

	if _polap_contains_step 5 "${_step_array[@]}"; then
		_polap_log0 "step 5"
		local mt_gfa="${_arg_long_reads%.*}.mt.gfa"
		cp -p "${_arg_outdir}/mt.1.gfa" "${mt_gfa}"
		_polap_log0 "output assembly graph: ${mt_gfa}"
	fi

	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
