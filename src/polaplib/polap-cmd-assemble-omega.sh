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

function _run_polap_assemble-omega {
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
  polap assemble-w - assemble after adjusting w option so that subsampling rate is between 0.1 and 0.5.

Synopsis:
  polap assemble-w [options]

Description:
  Rewrite this:

  polap assemble-rate does assemble organelle genomes by finding a sampling
  rate of input long-read data so that mapped reads are about 50x the seed
  contig size.
  Seed conitgs are very similar to a target genome; e.g., plastid genome
  assemble graph. When reads mapped on the seed contigs are over 50x the seed
  contig size, we sample reads so that reads are upto 50x times the seed contig
  size so that Flye genome assembly can be efficient for time and memory.
  However, 50x sampling reads can result in too small rate if mapped reads are
  way too much because of too much input long-read data. If the sampling rate of
  mapped reads is too small or less then 0.01, the randomness can lead to
  improper sampling, which results in worse genome assembly than the seed.
  In short, even though we assemble a good plastid-like graph for a seed,
  the second assembly from the seed can be a incomplete plastid graph because
  of too small sampling rate. Therefore, Input long-read data need to be
  reduced or subsampled before we map them on the seed contig so
  that sampling rate of mapped read is not too small.
  We start with 1 Gb of the input long-read data.
  Let p be a sampling rate of mapped reads. If p is less than 0.1, we reduce
  the input long-read by the rate of 0.1/p. If p is greater than 0.5,
  we increase it by two-fold if there are data available.
  We repeat this adjustemnt of subsampling rate of the input data so that
  p be between 0.1 and 0.5.

Options:
  -l FASTQ
    reads data file

  -i STR
    index of the source Flye assembly

  -j STR
    index of the destination Flye assembly

  -a FLOAT
    the lower bound of a subsampling rate range

  -b FLOAT
    the upper bound of a subsampling rate range

  -t INT
    the number of CPU cores

  -c INT
    maximum coverage of reads 

Examples:
  Get organelle genome sequences:
    polap assemble-w -l l.fq

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

	if ! _polap_gfatools-gfa2fasta; then
		_polap_error_message $?
	fi

	check_file_existence "${_polap_var_mtcontigname}"
	check_file_existence "${_polap_var_ga_contigger_edges_fasta}"

	if [[ ! -s "${_polap_var_mtcontigname}" ]]; then
		_polap_log0 "ERROR: no such file: ${_polap_var_mtcontigname}"
	fi

	_polap_lib_conda-ensure_conda_env polap || exit 1

	_run_polap_map-reads

	_polap_lib_oga-estimate-omega

	# remove NUMT/NUPT using rkmerrc
	# if [[ "${_arg_data_type}" == "pacbio-hifi" ]] && [[ "${_arg_plastid}" == "off" ]]; then
	if [[ "${_arg_data_type}" == "pacbio-hifi" ]]; then
		local _pread_sel="ptgaul-reads"
		local index=$(<"${_polap_var_oga_contig}/index.txt")
		local fq="${_polap_var_oga_seeds}/${_pread_sel}/${index}.fq"
		gunzip "${fq}.gz"
		local PREFIX="${_arg_outdir}/kmer/rmkc"
		local CLEANED="$PREFIX.cleaned.fastq.gz"
		_polap_log0 "rmkc on ${fq}"
		_polap_filter-reads-by-rmkc "${fq}"
		_polap_log0 "rmkc produces ${CLEANED}"
		_polap_lib_oga-flye-select-reads "${CLEANED}"
	else
		# echo _polap_lib_oga-flye-select-reads
		_arg_redo="on"
		_run_polap_assemble2
	fi

	conda deactivate

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
