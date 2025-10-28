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

function _run_polap_recruit {
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
what this main driver does (quick checklist)
	•	logs to its own recruit.log inside -o/--outdir.
	•	builds the plastid two-isomer panel first (polap-bash-pt-twoforms-from-IR-SSC.sh) so competitive mapping always has the “correct” chloroplast presentation to subtract.
	•	runs competitive mapping once (round0) against mito seeds ∪ plastid isomers ∪ optional nuclear decoy; only reads whose best hit is the mito target and whose aligned span ≥ TSPAN_MIN (preset-dependent) survive.
	•	optionally runs the k-mer signature recruiter and/or read-overlap recruiter, each with the 5% incremental stop rule pre-wired and controlled by --rounds.
	•	ends with a coverage sanity step if polap-bash-cov-quickcheck.sh is present (uses minimap2 -> samtools to plot/TSV summary).

preset knobs you can tweak
	•	conservative (default): higher identity, span, qcov — safest when drift is a concern.
	•	balanced: middle-ground defaults; good starting point for many datasets.
	•	aggressive-intergenic: lower id, qcov, lower TSPAN_MIN, increased MinHash scaled (i.e., sparser sketches) and shorter overlap minimum — this is meant to pull more noisy & intergenic reads [<0;162;19M[<0;162;19m[<0;155;18M[<0;155;18mwhile relying on competitive mapping plus span to hold off drift.

Options:
  -l FASTQ
    reads data file

Examples:
  Get organelle genome sequences:
    polap get-mtdna --plastid --species "Trifolium pratense"
    polap ${polap_cmd} 

  Get mitochondrial rRNA sequences from NCBI:
    polap ${polap_cmd} rrna

  Try assemble mtDNA:
    polap ${polap_cmd} --preset hifi

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (2024-2025)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_arg_menu[1]} == "help" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	if [[ "${_arg_help}" == "on" ]]; then
		bash "${_POLAPLIB_DIR}/polap-bash-recruit-mt-main.sh" --help >&3
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

	if [[ "${_arg_menu[1]}" == "busco" ]]; then
		mkdir -p refs
		_polap_lib_conda-ensure_conda_env polap-ncbitools || exit 1
		bash "${_POLAPLIB_DIR}/polap-bash-fetch-mt-rrna-and-busco.sh" -o refs -v
		conda deactivate
	fi

	if [[ "${_arg_menu[1]}" == "rrna" ]]; then
		_polap_log0 "downloading mitochondrial rRNA from NCBI ..."
		mkdir -p refs
		_polap_lib_conda-ensure_conda_env polap-ncbitools || exit 1
		bash "${_POLAPLIB_DIR}/polap-bash-fetch-mt-markers.sh" -o refs -v
		conda deactivate
		_polap_log0 "copy refs/plant_mt_rRNA.fna ${_POLAPLIB_DIR}/polap-mt.rrna.1.fna"
		_polap_log0 "file: refs/plant_mt_rRNA.fna"
	fi

	if [[ "${_arg_menu[1]}" == "infile" ]]; then
		_polap_lib_conda-ensure_conda_env polap || exit 1
		_polap_log0 bash "${_POLAPLIB_DIR}/polap-bash-recruit-mt-main.sh" \
			-i "${_arg_long_reads}" \
			-o "${_arg_outdir}/recruit" \
			--threads "${_arg_threads}" \
			--mt-seed "${_arg_outdir}/mt.0.fa" \
			--pt-ref "${_arg_unpolished_fasta}" \
			--seed-cds "${_POLAPLIB_DIR}/polap-mt.1.c70.3.fna" \
			--mt-rrna "${_POLAPLIB_DIR}/polap-mt.rrna.1.fna" \
			--target-bases 11g \
			--kmer-tool sourmash \
			-v
		conda deactivate
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
