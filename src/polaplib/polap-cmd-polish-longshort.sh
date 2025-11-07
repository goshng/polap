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

function _run_polap_polish-longshort {
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
  polap ${polap_cmd} - polish using Racon, fmlrc2, and Polypolish

Synopsis:
  polap ${polap_cmd} [options]

Description:
  polap ${polap_cmd} uses plastid and organelle genes to annotate reads
  using minimap2.

Options:
  -l FASTQ
    reads data file

  -a FASTQ
    short-read file

  -b FASTQ
    short-read file

  --infile1 FASTA
    unpolished mito genome sequence file

  --infile2 FASTA
    unpolished plastid genome sequence file

  --outfile FASTA
    polished genome sequence file

Examples:
  Polish mitochondrial genome sequences:
    polap ${polap_cmd} -l l.fq -a s1.fq -b s2.fq --infile1 mt.0.fa --infile2 pt.0.fa --outfile mt.1.fa

  Polish plastid genome sequences:
    polap ${polap_cmd} --plastid -l l.fq -a s1.fq -b s2.fq --infile1 mt.0.fa --infile2 pt.0.fa --outfile pt.1.fa

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

	_polap_lib_conda-ensure_conda_env polap-polish || exit 1

	if command -v racon >/dev/null 2>&1; then
		_polap_log1 "[Info] found racon, but not sure if it works."
	else
		_polap_log0 "[Info] racon is missing; the conda version of racon fails."
		_polap_log0 "[Info] execute: bolap setup-racon"
	fi

	local OUTDIR="${_arg_outdir}/polish-longshort"
	_arg_plastid="on"
	if [[ "${_arg_plastid}" == "on" ]]; then
		OUTDIR="${_arg_outdir}/polish-longshort/pt"
	else
		OUTDIR="${_arg_outdir}/polish-longshort/mt"
	fi
	_polap_log3_cmd rm -rf "$OUTDIR"

	if [[ "${_arg_long_reads_is}" == "on" && "${_arg_short_read1_is}" == "on" ]]; then

		# _polap_log3_cmd bash "${_POLAPLIB_DIR}/polap-bash-polish-hybrid-racon-fmlrc2-polypolish.sh" \
		# 	--ont "${_arg_long_reads}" \
		# 	--sr1 "${_arg_short_read1}" \
		# 	--sr2 "${_arg_short_read2}" \
		# 	--fasta "${_arg_infile}" \
		# 	--outdir "${OUTDIR}" \
		# 	--threads "${_arg_half_threads}" \
		# 	--rounds 2 --min-ident 0.80 --min-alen 2000 \
		# 	${_arg_verbose_str}

		if [[ "${_arg_plastid}" == "on" ]]; then
			_polap_log3_cmd bash "${_POLAPLIB_DIR}/polap-bash-polish-racon-polypolish.sh" \
				--target cp \
				--fasta "${_arg_infile1}" \
				--other "${_arg_infile2}" \
				--ont "${_arg_long_reads}" \
				--sr1 "${_arg_short_read1}" \
				--sr2 "${_arg_short_read2}" \
				--outdir "${OUTDIR}" \
				--threads "${_arg_half_threads}" \
				--racon-rounds 3
		else
			# Paste the scripts above into the touched paths, then run:
			_polap_log3_cmd bash "${_POLAPLIB_DIR}/polap-bash-polish-racon-polypolish.sh" \
				--target mt \
				--fasta "${_arg_infile1}" \
				--other "${_arg_infile2}" \
				--ont "${_arg_long_reads}" \
				--sr1 "${_arg_short_read1}" \
				--sr2 "${_arg_short_read2}" \
				--outdir "${OUTDIR}" \
				--threads "${_arg_half_threads}" \
				--racon-rounds 3 \
				--mask auto --mask-min-ident 0.85 --mask-min-len 150 --mask-pad 500 --mask-qc 1 \
				--cov-bin 200 --cov-min-mapq 0

		fi

	elif [[ "${_arg_long_reads_is}" == "on" ]]; then

		bash "${_POLAPLIB_DIR}/polap-bash-polish-hybrid-racon-fmlrc2-polypolish.sh" \
			--ont "${_arg_long_reads}" \
			--fasta "${_arg_infile}" \
			--outdir "${OUTDIR}" \
			--threads "${_arg_half_threads}" \
			--rounds 2 --min-ident 0.80 --min-alen 2000 \
			${_arg_verbose_str}

	elif [[ "${_arg_short_read1_is}" == "on" ]]; then

		bash "${_POLAPLIB_DIR}/polap-bash-polish-hybrid-racon-fmlrc2-polypolish.sh" \
			--sr1 "${_arg_short_read1}" \
			--sr2 "${_arg_short_read2}" \
			--fasta "${_arg_infile}" \
			--outdir "${OUTDIR}" \
			--threads "${_arg_half_threads}" \
			--rounds 2 --min-ident 0.80 --min-alen 2000 \
			${_arg_verbose_str}

	else
		_polap_log0 "[ERROR] Reads are required."
	fi

	_polap_log1_cmd cp -p "${OUTDIR}/polished.fa" "${_arg_outfile}"

	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
