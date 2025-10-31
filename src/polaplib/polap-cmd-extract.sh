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

function _run_polap_extract {
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

  --infile1
    plastid gfa

  --infile2
    mitochondrial gfa

Examples:
  Get organelle genome sequences:
    polap ${polap_cmd} --infile1 pt.1.gfa --infile2 mt.1.gfa

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

	_polap_lib_conda-ensure_conda_env polap-oatk || exit 1

	local OUTDIR="${_arg_outdir}/extract"
	_polap_log3_cmd rm -rf "$OUTDIR"
	_polap_log3_cmd mkdir -p "$OUTDIR"

	_polap_log1 "Normalize PT Flye's assembly gfa file"
	bash "${_POLAPLIB_DIR}/polap-bash-gfa-normalize.sh" \
		--gfa "${_arg_infile1}" \
		--out "${OUTDIR}/pt.norm.gfa" \
		--idmap "${OUTDIR}/pt.idmap.tsv" \
		--ec-mode mean --ec-round round 2>"${_polap_output_dest}"

	_polap_log1 "Normalize MT Flye's assembly gfa file"
	# _polap_log0 bash "${_POLAPLIB_DIR}/polap-bash-gfa-normalize.sh" \
	# 	--gfa "${_arg_infile2}" \
	# 	--out "${OUTDIR}/mt.norm.gfa" \
	# 	--idmap "${OUTDIR}/mt.idmap.tsv" \
	# 	--ec-mode mean --ec-round round
	bash "${_POLAPLIB_DIR}/polap-bash-gfa-normalize.sh" \
		--gfa "${_arg_infile2}" \
		--out "${OUTDIR}/mt.norm.gfa" \
		--idmap "${OUTDIR}/mt.idmap.tsv" \
		--ec-mode mean --ec-round round 2>"${_polap_output_dest}"

	_polap_log1 "Oatk's annotation on PT assembly"
	local OATK_HMMANNOT="hmmannot"
	if command -v hmmannot >/dev/null 2>&1; then
		OATK_HMMANNOT="hmmannot"
	elif command -v hmm_annotation >/dev/null 2>&1; then
		OATK_HMMANNOT="hmm_annotation"
	else
		_polap_log0 "[ERROR] neither hmmannot nor hmm_annotation found."
		exit 127
	fi

	$OATK_HMMANNOT -t 8 -o "${OUTDIR}/annot_pltd.txt" \
		OatkDB/v20230921/embryophyta_pltd.fam \
		"${OUTDIR}/pt.norm.gfa" 2>"${_polap_output_dest}"

	_polap_log1 "Oatk's annotation on MT assembly"
	$OATK_HMMANNOT -t 8 -o "${OUTDIR}/annot_mito.txt" \
		OatkDB/v20230921/embryophyta_mito.fam \
		"${OUTDIR}/mt.norm.gfa" 2>"${_polap_output_dest}"

	_polap_log1 "Oatk's pathfinder in PT assembly for ptDNA sequence"
	pathfinder -p "${OUTDIR}/annot_pltd.txt" \
		-o "${OUTDIR}/oatk" \
		"${OUTDIR}/pt.norm.gfa" 2>"${_polap_output_dest}"

	_polap_log1 "Oatk's pathfinder in MT assembly for ptDNA sequence"
	pathfinder -m "${OUTDIR}/annot_mito.txt" \
		-o "${OUTDIR}/oatk" \
		"${OUTDIR}/mt.norm.gfa" 2>"${_polap_output_dest}"

	_polap_log0 "ptDNA: ${OUTDIR}/oatk.pltd.ctg.fasta"
	_polap_log0 "mtDNA: ${OUTDIR}/oatk.mito.ctg.fasta"

	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
