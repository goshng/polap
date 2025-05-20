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
# Tip!
# How to extract commands that were expected:
# src/polap.sh reduce-data --redo -v -v -v 2>&1 | grep -E "^rm|^seqkit|^ln"
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

source "${_POLAPLIB_DIR}/polap-function-set-variables.sh"
source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"

################################################################################
# Runs the whole-genome assembly.
################################################################################

function _polap_summary-generated-reads {
	local _fastq_gz=$1
	local _fastq_stats=$2

	if [ -s "${_fastq_stats}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log1 "  found: ${_fastq_stats}, so skipping POLAP long-read statisics ..."
		_polap_log2_cat "${_fastq_stats}"
	else
		if [ -s "${_fastq_gz}" ]; then
			_polap_log3_pipe "seqkit stats -T ${_fastq_gz} \
				>${_fastq_stats}"
			_polap_log1 "  output: ${_fastq_stats}"
			_polap_log2_cat "${_fastq_stats}"
		else
			_polap_log1 "  no such file: ${_fastq_gz}, so skipping read statisics ..."
		fi
	fi
}

################################################################################
# Statisics of the short-read and POLAP's long-read (nk.fq.gz) dataset.
################################################################################
function _run_polap_summary-reads { # statisics of the read dataset
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set paths for bioproject data
	_polap_set-variables-short-read
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
# Summarize the statistics of the long- and short-read dataset.
#
# Arguments:
#   -a ${_arg_short_read1}: a short-read fastq data file
#   -b ${_arg_short_read2}: another short-read fastq data file (optional: not tested yet!)
#   or
#   --bioproject use
#   -o ${_arg_outdir}
# Inputs:
#   ${_arg_short_read1}: a short-read fastq data file
#   ${_arg_short_read2}: another short-read fastq data file (optional: not tested yet!)
#   ${_polap_var_outdir_nk_fq_gz}: POLAP generate nk.fq.gz
#   ${_polap_var_outdir_lk_fq_gz}: POLAP generate lk.fq.gz
# Outputs:
#   ${_polap_var_outdir_s1_fq_stats}: short-read data statisics
#   ${_polap_var_outdir_s2_fq_stats}: short-read data statisics
#   ${_polap_var_outdir_nk_fq_stats}: POLAP long-read data statisics
#   ${_polap_var_outdir_lk_fq_stats}: POLAP long-read data statisics
# Precondition:
#   (for BioProjectID case)
#   get-bioproject --bioproject <BioProjectID> -o ${_arg_outdir}
Example: $(basename "$0") ${_arg_menu[0]} -a <file> [-b <file>]
Example: $(basename "$0") ${_arg_menu[0]} -o ${_arg_outdir} --bioproject use
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [ -s "${_polap_var_outdir_s1_fq_stats}" ]; then
			_polap_log0_cat "${_polap_var_outdir_s1_fq_stats}"
		else
			_polap_log0 "No short-read s1 statisics"
		fi

		if [ -s "${_polap_var_outdir_s2_fq_stats}" ]; then
			_polap_log0_cat "${_polap_var_outdir_s2_fq_stats}"
		else
			_polap_log0 "No short-read s2 statisics"
		fi

		if [ -s "${_polap_var_outdir_nk_fq_stats}" ]; then
			_polap_log0_cat "${_polap_var_outdir_nk_fq_stats}"
		else
			_polap_log0 "No POLAP long-read statisics"
		fi

		if [ -s "${_polap_var_outdir_lk_fq_stats}" ]; then
			_polap_log0_cat "${_polap_var_outdir_lk_fq_stats}"
		else
			_polap_log0 "No POLAP long-read statisics"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	_polap_log0 "computing stats of the short-read sequencing data ..."
	if [[ -d "${_arg_outdir}" ]]; then
		_polap_log2 "  output folder: ${_arg_outdir}"
	else
		_polap_log3_cmd mkdir -p "${_arg_outdir}"
	fi

	# check_file_existence "${_arg_short_read1}"
	# _polap_log1 "  input1: ${_arg_short_read1}"
	# if [ -s "${_arg_short_read2}" ]; then
	# 	_polap_log1 "  input2: ${_arg_short_read2}"
	# else
	# 	_polap_log1 "  input2: no such file: ${_arg_short_read2}"
	# 	_polap_log1 "    we use a single short-read data file."
	# fi
	#
	# if [ -s "${_polap_var_outdir_fq_stats}" ] && [ "${_arg_redo}" = "off" ]; then
	# 	_polap_log1 "  found: ${_polap_var_outdir_fq_stats}, so skipping short-read statisics ..."
	# else
	# 	if [[ -s "${_arg_short_read1}" ]]; then
	#
	# 		_polap_log3_pipe "seqkit stats -T ${_arg_short_read1} |\
	#      csvtk del-header >${_polap_var_outdir_fq_stats}"
	#
	# 		if [ -s "${_arg_short_read2}" ]; then
	# 			_polap_log3_pipe "seqkit stats -T ${_arg_short_read2} |\
	# 			csvtk del-header >>${_polap_var_outdir_fq_stats}"
	# 		fi
	# 	else
	# 		_polap_log0 "  no such file: ${_arg_short_read1}"
	# 	fi
	# fi
	#
	# _polap_log1 "  output1: ${_polap_var_outdir_fq_stats}"
	# _polap_log2_cat "${_polap_var_outdir_fq_stats}"

	_polap_summary-generated-reads \
		"${_arg_long_reads}" \
		"${_polap_var_outdir_l_fq_stats}"

	_polap_summary-generated-reads \
		"${_polap_var_outdir_nk_fq_gz}" \
		"${_polap_var_outdir_nk_fq_stats}"

	_polap_summary-generated-reads \
		"${_polap_var_outdir_lk_fq_gz}" \
		"${_polap_var_outdir_lk_fq_stats}"

	# Get the total long-read size
	local _l_size=$(awk 'NR==2 { print $5 }' "${_polap_var_outdir_l_fq_stats}")
	local _lk_size=$(awk 'NR==2 { print $5 }' "${_polap_var_outdir_lk_fq_stats}")
	local _nk_size=$(awk 'NR==2 { print $5 }' "${_polap_var_outdir_nk_fq_stats}")
	local _EXPECTED_GENOME_SIZE=$(<"${_polap_var_outdir_genome_size}")
	_EXPECTED_GENOME_SIZE=${_EXPECTED_GENOME_SIZE%.*}
	local _l_coverage=$(echo "scale=3; ${_l_size}/${_EXPECTED_GENOME_SIZE}" | bc)
	local _lk_coverage=$(echo "scale=3; ${_lk_size}/${_EXPECTED_GENOME_SIZE}" | bc)
	local _nk_coverage=$(echo "scale=3; ${_nk_size}/${_EXPECTED_GENOME_SIZE}" | bc)

	{
		echo "genome size: ${_EXPECTED_GENOME_SIZE}"
		echo "long-read size: ${_l_size}"
		echo "long-read coverage: ${_l_coverage}"
		echo "lk size: ${_lk_size}"
		echo "lk coverage: ${_lk_coverage}"
		echo "nk size: ${_nk_size}"
		echo "nk coverage: ${_nk_coverage}"
	} >"${_polap_var_outdir_fq_stats}"

	_polap_log1 NEXT: $(basename "$0") total-length-long -o "${_arg_outdir}" -l ${_arg_long_reads}
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Computes the total number of nucleotides of long-read data.
#
# Arguments:
#   -l $_arg_long_reads: a long-read fastq data file
# Inputs:
#   $_arg_long_reads: a long-read fastq data file
#   ${_arg_outdir}
# Outputs:
#   ${_arg_outdir}/long_total_length.txt
################################################################################
function _run_polap_total-length-long { # total size (bp) of long-read data
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# CHECK: local function
	_polap_set-variables-long-read
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"
	source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"

	help_message=$(
		cat <<HEREDOC
Compute the total number of nucleotides of long-read data.

Inputs
------

- long-read fastq data file: ${_arg_long_reads}

Outputs
-------

- ${_polap_var_outdir_long_total_length}

Arguments
---------

-o ${_arg_outdir}
-l ${_arg_long_reads}: a long-read fastq data file

# Arguments:
#   -l ${_arg_long_reads}: a long-read fastq data file (the highest priority)
#   or
#   --bioproject use
#   -o ${_arg_outdir}: ${_arg_outdir}/0-bioproject (the least priority)
#

Precondition
------------
(for BioProjectID case)
get-bioproject --bioproject <BioProjectID> -o ${_arg_outdir}

Usages
------
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads}
$(basename "$0") ${_arg_menu[0]} -o ${_arg_outdir} --bioproject use
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [ -s "${_polap_var_outdir_long_total_length}" ]; then
			_polap_log0_cat "${_polap_var_outdir_long_total_length}"
		else
			_polap_log0 "No long-read total length"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	_polap_log1 "  count the total number of bases in the long-read dataset of $_arg_long_reads ..."
	if [[ ! -d "${_arg_outdir}" ]]; then
		_polap_log2 "  creating output folder: ${_arg_outdir}"
		_polap_log3_cmd mkdir -p "${_arg_outdir}"
	fi
	check_file_existence "${_arg_long_reads}"

	_polap_log2 "  input: ${_arg_long_reads}"

	if [ -s "${_polap_var_outdir_long_total_length}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  found: ${_polap_var_outdir_long_total_length}, so skipping long-read statisics ..."
	else
		_polap_log3_pipe "seqkit stats -Ta ${_arg_long_reads} |\
			csvtk cut -t -f sum_len |\
			csvtk del-header \
				>${_polap_var_outdir_long_total_length}"
	fi

	_polap_log2 "  output: ${_polap_var_outdir_long_total_length}"
	_polap_log3_cat "${_polap_var_outdir_long_total_length}"

	local LONG_TOTAL_LENGTH=$(<"${_polap_var_outdir_long_total_length}")
	local _total_long_read=$(_polap_utility_convert_bp ${LONG_TOTAL_LENGTH})
	_polap_log1 "  total length of the long-read dataset (bases): ${_total_long_read}"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_total-length-short { # total size (bp) of short-read data
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# CHECK: local function
	_polap_set-variables-long-read
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"
	source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"

	help_message=$(
		cat <<HEREDOC
Compute the total number of nucleotides of short-read data.

Inputs
------

- short-read fastq data file: ${_arg_short_read1}
- short-read fastq data file: ${_arg_short_read2}

Outputs
-------

- ${_polap_var_outdir_short_total_length}

Arguments
---------

-o ${_arg_outdir}
-a ${_arg_short_read1}: a short-read fastq data file

Usages
------
$(basename "$0") ${_arg_menu[0]} -a ${_arg_short_read1}
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [ -s "${_polap_var_outdir_short_total_length}" ]; then
			_polap_log0_cat "${_polap_var_outdir_short_total_length}"
		else
			_polap_log0 "No short-read total length"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	_polap_log0 "counting the total number of bases in the short-read dataset ..."
	if [[ ! -d "${_arg_outdir}" ]]; then
		_polap_log2 "  creating output folder: ${_arg_outdir}"
		_polap_log3_cmd mkdir -p "${_arg_outdir}"
	fi

	check_file_existence "${_arg_short_read1}"

	_polap_log1 "  input: ${_arg_short_read1}"

	if [ -s "${_polap_var_outdir_short_total_length}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log1 "  found: ${_polap_var_outdir_short_total_length}, so skipping short-read statisics ..."
	else
		_polap_log3_pipe "seqkit stats -Ta ${_arg_short_read1} |\
			csvtk cut -t -f sum_len |\
			csvtk del-header \
				>${_polap_var_outdir_short_total_length}"
	fi

	_polap_log1 "  output: ${_polap_var_outdir_short_total_length}"
	_polap_log2_cat "${_polap_var_outdir_short_total_length}"

	local _TOTAL_LENGTH=$(<"${_polap_var_outdir_short_total_length}")
	local _total_short_read=$(_polap_utility_convert_bp ${_TOTAL_LENGTH})
	_polap_log0 "  total length of the short-read dataset (bases): ${_total_short_read}"

	_polap_log1 NEXT menu: find-genome-size

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Estimates the whole genome size using short-read data.
#
# Steps:
#   1. jellyfish count
#   2. jellyfish histo
#   3. estimate the genome siz.
#
# Arguments:
#   -a ${_arg_short_read1}: a short-read fastq data file
#   -b ${_arg_short_read2}: another short-read fastq data file (optional)
#
# Outputs:
#   ${_arg_outdir}/jellyfish_out
#   ${_arg_outdir}/jellyfish_out.histo
#   ${_arg_outdir}/short_expected_genome_size.txt
################################################################################
function _run_polap_find-genome-size { # estimate the whole genome size
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set paths for bioproject data
	_polap_set-variables-short-read
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'
	source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"

	help_message=$(
		cat <<HEREDOC
Estimate the whole genome size applying the short-read dataset to JellyFish.
#
# Arguments:
#   -a ${_arg_short_read1}: a short-read fastq data file
#   -b ${_arg_short_read2}: another short-read fastq data file (optional)
#   or
#   --bioproject use
#   -o ${_arg_outdir}
# Inputs:
#   ${_arg_short_read1}: a short-read fastq data file
#   ${_arg_short_read2}: another short-read fastq data file (optional)
# Outputs:
#   ${_polap_var_outdir_genome_size}
Example: $(basename "$0") ${_arg_menu[0]} -a <file> [-b <file>]
Example: $(basename "$0") ${_arg_menu[0]} -o ${_arg_outdir} --bioproject use
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [ -s "${_polap_var_outdir_genome_size}" ]; then
			_polap_log0_file "${_polap_var_outdir_genome_size}"
			_polap_log0_cat "${_polap_var_outdir_genome_size}"
		else
			_polap_log0 "No genome size estimate from a short-read dataset"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "estimating the genome size using your short-read data [${_arg_short_read1}] and [${_arg_short_read2}] ..."

	if [[ -d "${_arg_outdir}" ]]; then
		_polap_log2 "  output folder: ${_arg_outdir}"
	else
		_polap_log3 mkdir -p "${_arg_outdir}"
		mkdir -p "${_arg_outdir}"
	fi
	check_file_existence "${_arg_short_read1}"

	_polap_log1 "  input1: ${_arg_short_read1}"
	if [ -s "${_arg_short_read2}" ]; then
		_polap_log1 "  input2: ${_arg_short_read2}"
	else
		_polap_log1 "  input2: no such file: ${_arg_short_read2}"
		_polap_log1 "    we use a single short-read data file."
	fi

	# See https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/
	if [ -s "${_polap_var_outdir_genome_size}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log1 "  found: ${_polap_var_outdir_genome_size}, so skipping the genome size estimation using the short-read data ..."
		_polap_log1_file "${_polap_var_outdir_genome_size}"
	else
		if [ -s "${_polap_var_outdir_jellyfish_out}" ] && [ "${_arg_redo}" = "off" ]; then
			_polap_log1 "  found: ${_polap_var_outdir_jellyfish_out}, so skipping the JellyFish count step using the short-read data ..."
			_polap_log1_file "${_polap_var_outdir_jellyfish_out}"
		else
			if [[ -s "${_arg_short_read1}" ]] && [[ -s "${_arg_short_read2}" ]]; then
				_polap_log3_pipe "jellyfish count \
					-t ${_arg_threads} -C -m 19 \
					-s 5G \
					-o ${_polap_var_outdir_jellyfish_out} \
					--min-qual-char=? \
					${_arg_short_read1} ${_arg_short_read2}"
			elif [[ -s "${_arg_short_read1}" ]] && [[ ! -s "${_arg_short_read2}" ]]; then
				_polap_log3_pipe "jellyfish count \
					-t ${_arg_threads} -C -m 19 \
					-s 5G \
					-o ${_polap_var_outdir_jellyfish_out} \
					--min-qual-char=? \
					${_arg_short_read1}"
			else
				die "ASSERT: we must have at least one short-read fastq file."
			fi
			check_file_existence "${_polap_var_outdir_jellyfish_out}"
		fi

		check_file_existence "${_polap_var_outdir_jellyfish_out}"
		_polap_log3_pipe "jellyfish histo \
			-o ${_polap_var_outdir_jellyfish_out_histo} \
			${_polap_var_outdir_jellyfish_out}"
		_polap_log2_file "${_polap_var_outdir_jellyfish_out_histo}"

		_polap_log3_pipe "Rscript ${_POLAPLIB_DIR}/run-polap-r-jellyfish.R \
			${_polap_var_outdir_jellyfish_out_histo} \
			${_polap_var_outdir_genome_size}"
		# Check the exit status
		if [ $? -ne 0 ]; then
			# Take action if needed, e.g., logging, sending a notification, etc.
			die "ERROR: not being to estimate the genome size using short-read data: short-read may be too small."
		fi
	fi
	_polap_log1 "  output: ${_polap_var_outdir_genome_size}"
	_polap_log2_cat "${_polap_var_outdir_genome_size}"

	local _EXPECTED_GENOME_SIZE=$(<"${_polap_var_outdir_genome_size}")
	local _EXPECTED_GENOME_SIZE=${_EXPECTED_GENOME_SIZE%.*}
	local _expected_genome_size_bp=$(_polap_utility_convert_bp ${_EXPECTED_GENOME_SIZE})
	_polap_log0 "  expected genome size using short-read data (bases): ${_expected_genome_size_bp}"

	_polap_log0 "  summary statisics of the short-read data ..."
	if [[ -s "${_polap_var_outdir_s1_fq_stats}" ]]; then
		_polap_log1 "  found: ${_polap_var_outdir_s1_fq_stats}"
	else
		_polap_summary-generated-reads \
			"${_arg_short_read1}" \
			"${_polap_var_outdir_s1_fq_stats}"
	fi

	if [[ -s "${_polap_var_outdir_s2_fq_stats}" ]]; then
		_polap_log1 "  found: ${_polap_var_outdir_s2_fq_stats}"
	else
		_polap_summary-generated-reads \
			"${_arg_short_read2}" \
			"${_polap_var_outdir_s2_fq_stats}"
	fi

	_polap_log1 NEXT: $(basename "$0") reduce-data -o "${_arg_outdir}" -l "${_arg_long_reads}" [-m "${_arg_min_read_length}"]

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Steps:
# 1. Firstly, generate a long-read DNA sequence of ${_arg_min_read_length}
#    base pairs (bp) to serve as the foundation for further analysis and
#    experimentation.
# 2. subsample the long-read data to achieve a desired level of coverage.
#    Checks whether the long-read coverage falls below the specified threshold of ${_arg_coverage}.
#    If the condition is met, retain the full-length read data;
#    otherwise, sample the long reads up to a certain level of coverage.
# 3. delete long reads that are shorter than a specified sequence length
#    threshold, such as 3 kilobases.
################################################################################
function _run_polap_reduce-data { # reduce the long-read data, if too big
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# CHECK: local function
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"
	_polap_set-variables-long-read
	source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"

	help_message=$(
		cat <<HEREDOC
This tool reduces long-read data for a Flye genome assembly.

Inputs
------

- source assembly number: i
- ${_polap_var_outdir_genome_size}
- ${_polap_var_outdir_long_total_length}
- ${_arg_long_reads}

- destination assembly number: j
- j/01-contig/contig1.fa
- j/01-contig/contig2.fa
- ${_polap_var_outdir_lk_fq_gz}

Outputs
-------

- ${_polap_var_outdir_lk_fq_gz}: reads longer than ${_arg_min_read_length} bp
- ${_polap_var_outdir_nk_fq_gz}: subsample of ${_arg_long_reads}
- ${_polap_var_outdir_nk8_fq_gz}: subsample of ${_arg_long_reads}

Arguments
---------
-l ${_arg_long_reads}: a long-read fastq data file
or
--bioproject use
-o ${_arg_outdir}: ${_arg_outdir}/0-bioproject (the least priority)
-m ${_arg_min_read_length}: the long-read sequence length threshold
-c ${_arg_coverage}: the target coverage
--no-reduction-reads
--random-seed 11: for seqkit default seed
--test

Menu
----
split <N>: splits ${_polap_var_outdir_lk_fq_gz} into N parts for faster read selections

Usage
-----
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads}
$(basename "$0") ${_arg_menu[0]} -o ${_arg_outdir} --bioproject use
$(basename "$0") ${_arg_menu[0]} split 10
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [ -s "${_polap_var_outdir_nk_fq_gz}" ]; then
			_polap_log0_cat "${_polap_var_outdir_nk_fq_stats}"
		else
			_polap_log0 "No reduced subsample long-read data file: ${_polap_var_outdir_nk_fq_gz}"
		fi

		if [ -s "${_polap_var_outdir_lk_fq_gz}" ]; then
			_polap_log0_cat "${_polap_var_outdir_lk_fq_stats}"
		else
			_polap_log0 "No reduced long-read: ${_polap_var_outdir_lk_fq_gz}"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	if [[ "${_arg_menu[1]}" == "split" ]]; then
		if [ -s "${_polap_var_outdir_lk_fq_gz}" ]; then
			seqkit split2 -p "${_arg_menu[2]}" "${_polap_var_outdir_lk_fq_gz}"
		else
			_polap_log0 "No reduced long-read: ${_polap_var_outdir_lk_fq_gz}"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "reducing the long-read data ${_arg_long_reads} for the whole- and organelle-genome assemblies ..."

	# Check for required files
	check_folder_existence "${_arg_outdir}"
	check_file_existence "${_polap_var_outdir_long_total_length}"
	check_file_existence "${_polap_var_outdir_genome_size}"
	check_file_existence "${_arg_long_reads}"

	_polap_log1 "  input1: ${_arg_long_reads}"
	_polap_log1 "  input2: ${_polap_var_outdir_long_total_length}"
	_polap_log1 "  input3: ${_polap_var_outdir_genome_size}"

	# Get the expected genome size and long-read sequencing coverage
	local _LONG_TOTAL_LENGTH=$(<"${_polap_var_outdir_long_total_length}")
	local _long_total_length_bp=$(_polap_utility_convert_bp ${_LONG_TOTAL_LENGTH})
	local _EXPECTED_GENOME_SIZE=$(<"${_polap_var_outdir_genome_size}")
	local _expected_genome_size_bp=$(_polap_utility_convert_bp ${_EXPECTED_GENOME_SIZE})
	local _EXPECTED_LONG_COVERAGE=$(echo "scale=3; ${_LONG_TOTAL_LENGTH}/${_EXPECTED_GENOME_SIZE}" | bc)
	_polap_log2 "    LONG_TOTAL_LENGTH: ${_long_total_length_bp}"
	_polap_log2 "    EXPECTED_GENOME_SIZE: ${_expected_genome_size_bp}"
	_polap_log2 "    EXPECTED_LONG_COVERAGE: ${_EXPECTED_LONG_COVERAGE}x"

	if [[ -s "${_polap_var_outdir_nk_fq_gz}" ]] &&
		[[ -s "${_polap_var_outdir_lk_fq_gz}" ]] &&
		[[ "${_arg_redo}" == "off" ]]; then
		# [[ -s "${_polap_var_outdir_nk4_fq_gz}" ]] &&
		_polap_log0 "  found: ${_polap_var_outdir_nk_fq_gz}, so skipping the long-read data reduction."
		_polap_log0 "  found: ${_polap_var_outdir_lk_fq_gz}, so skipping the long-read data reduction."
		# _polap_log0 "  found: ${_polap_var_outdir_nk4_fq_gz}, so skipping the long-read data reduction."
		return
	fi

	_polap_log0 "  reducing ${_arg_long_reads} -> ${_polap_var_outdir_lk_fq_gz} with the minimum length of ${_arg_min_read_length} bp"

	if [[ -s "${_polap_var_outdir_lk_fq_gz}" ]] && [[ "${_arg_redo}" == "off" ]]; then
		_polap_log0 "  found: ${_polap_var_outdir_lk_fq_gz}, so skipping the long-read data reduction."
	else
		_polap_log1 "  keeps long reads of length being at least ${_arg_min_read_length} bp ..."
		_polap_log2 "    input1: ${_arg_long_reads}"
		_polap_log2 "    output: ${_polap_var_outdir_lk_fq_gz}"
		_polap_log3_cmd rm -f "${_polap_var_outdir_lk_fq_gz}"
		_polap_log3_pipe "seqkit seq \
      --quiet \
      -m ${_arg_min_read_length} \
      --threads 4 \
		  ${_arg_long_reads} \
		  -o ${_polap_var_outdir_lk_fq_gz} \
      >${_polap_output_dest} 2>&1"

		_polap_log1 "  creating the statisics for the size reduced long-read data ..."
		_polap_log2 "    input1: ${_polap_var_outdir_lk_fq_gz}"
		_polap_log2 "    output: ${_polap_var_outdir_lk_fq_stats}"
		_polap_log3_pipe "seqkit stats -aT \
      ${_polap_var_outdir_lk_fq_gz} \
      >${_polap_var_outdir_lk_fq_stats}"
		_polap_log2_column "${_polap_var_outdir_lk_fq_stats}"
	fi

	_polap_log0 "  reducing ${_arg_long_reads} -> ${_polap_var_outdir_nk_fq_gz} with the target coverage of ${_arg_coverage}x ..."

	if [[ -s "${_polap_var_outdir_nk_fq_gz}" ]] && [[ "${_arg_redo}" == "off" ]]; then
		_polap_log0 "  found: ${_polap_var_outdir_nk_fq_gz}, so skipping the long-read data reduction."
	else
		# subsample the long-read data so that the target coverage is ${_arg_coverage}.
		local nfq_file="${_arg_outdir}/n.fq"
		_polap_log2 "  deletes ${nfq_file} if there is one."
		_polap_log3_cmd rm -f "${nfq_file}"
		if [[ ${_arg_test} == "on" ]]; then
			_polap_log1 "  OPTION: --test : No reduction of the test long-read data"
			_polap_log3_cmd ln -s $(realpath "$_arg_long_reads") "$nfq_file"
		elif [[ ${_arg_reduction_reads} == "off" ]]; then
			_polap_log1 "  OPTION: --no-reduction-reads : No reduction of the long-read data"
			_polap_log3_cmd ln -s $(realpath "$_arg_long_reads") "$nfq_file"
		else
			# if [ "$EXPECTED_LONG_COVERAGE " -lt ${_arg_coverage} ]; then
			if echo "${_EXPECTED_LONG_COVERAGE} < ${_arg_coverage}" | bc -l | grep -q 1; then
				_polap_log1 "  No reduction of the long-read data because ${_EXPECTED_LONG_COVERAGE} < ${_arg_coverage}"
				_polap_log3_cmd ln -s $(realpath "$_arg_long_reads") "$nfq_file"
			else
				_polap_log1 "SUGGESTION: you might want to increase the minimum read lengths because you have enough long-read data."
				local _RATE=$(echo "scale=9; ${_arg_coverage}/${_EXPECTED_LONG_COVERAGE}" | bc)
				# Compare value with 0
				if echo "${_RATE} > 0" | bc -l | grep -q 1; then
					_polap_log1 "  sampling long-read data by ${_RATE} ..."
					_polap_log1 "    ${_RATE} <= target long-read genome coverage[${_arg_coverage}]/expected long-read genome coverage[${_EXPECTED_LONG_COVERAGE}] ..."
					_polap_lib_random-get
					local _random_seed=${_polap_var_random_number}
					_polap_log1 "  random seed for reducing the whole-genome assembly long-read data: ${_random_seed}"
					# _polap_log3 "seqkit sample -p ${_RATE} ${_arg_long_reads} -o ${nfq_file}"
					# seqkit sample -p "${_RATE}" "${_arg_long_reads}" -o "${nfq_file}" >${_polap_output_dest} 2>&1
					_polap_log3_pipe "seqkit sample \
            -p ${_RATE} \
            -s ${_random_seed} \
            ${_arg_long_reads} \
            2>${_polap_output_dest} |
            seqkit seq \
              --quiet \
              -m ${_arg_min_read_length} \
              --threads 4 \
		          -o ${_polap_var_outdir_nk_fq_gz} \
              >${_polap_output_dest} 2>&1"
					_polap_log3_pipe "echo ${_random_seed} >${_polap_var_outdir_nk_fq_gz}.random.seed.${_random_seed}"
					_polap_log1 "  ${nfq_file}: a reduced long-read data is created"
				else
					_polap_log0 "  target coverage: ${_arg_coverage}"
					_polap_log0 "  long-read coverage: ${_EXPECTED_LONG_COVERAGE}"
					_polap_log0 "  sampling rate is ${_arg_coverage} / ${_EXPECTED_LONG_COVERAGE} => ${_RATE}"
					_polap_log0 "  genome size: ${_EXPECTED_GENOME_SIZE}"
					_polap_log0 "  total long-read: ${_LONG_TOTAL_LENGTH}"
					_polap_log0 "  Too large expected long-read coverage"
					_polap_log0 "  Expected genome size may be too small."
					die "ERROR: long-read sampling rate is not greater than 0."
				fi
			fi
		fi

		if [[ -s "${nfq_file}" ]]; then
			# purge the long-read data of shorter than ${_arg_min_read_length} bp
			_polap_log1 "  keeps long reads of length being at least ${_arg_min_read_length} bp ..."
			_polap_log2 "    input1: ${nfq_file}"
			_polap_log2 "    output: ${_polap_var_outdir_nk_fq_gz}"
			_polap_log3_cmd rm -f "${_polap_var_outdir_nk_fq_gz}"
			_polap_log3_pipe "seqkit seq \
        --quiet \
        -m ${_arg_min_read_length} \
        --threads 4 \
		    ${nfq_file} \
  		  -o ${_polap_var_outdir_nk_fq_gz} \
        >${_polap_output_dest} 2>&1"

			_polap_log3_cmd rm -f "$nfq_file"
		fi

		_polap_log1 "  creating the statisics for the reduced long-read data ..."
		_polap_log3_pipe "seqkit stats -aT \
      ${_polap_var_outdir_nk_fq_gz} \
      >${_polap_var_outdir_nk_fq_stats}"
		_polap_log1 "    output: ${_polap_var_outdir_nk_fq_stats}"
		_polap_log2_column "${_polap_var_outdir_nk_fq_stats}"
	fi

	# nk4, nk8, nk16, nk32
	# _polap_log0 "  reducing ${_polap_var_outdir_lk_fq_gz} -> ${_polap_var_outdir_nk4_fq_gz} ..."
	#
	# if [[ -s "${_polap_var_outdir_nk4_fq_gz}" ]] && [[ "${_arg_redo}" == "off" ]]; then
	# 	_polap_log0 "  found: ${_polap_var_outdir_nk4_fq_gz}, so skipping the long-read data reduction."
	# else
	# 	local target_file=$(readlink -f "${_polap_var_outdir_lk_fq_gz}")
	# 	local file_size=$(stat --format="%s" "$target_file")
	# 	local limit_file_size=$((4 * 1024 * 1024 * 1024))
	# 	local ratio=$(echo "scale=4; $limit_file_size / $file_size" | bc)
	#
	# 	# if [ "$EXPECTED_LONG_COVERAGE " -lt ${_arg_coverage} ]; then
	# 	if echo "${file_size} < ${limit_file_size}" | bc -l | grep -q 1; then
	# 		_polap_log1 "  No reduction of the long-read data because ${_EXPECTED_LONG_COVERAGE} < ${_arg_coverage}"
	# 		_polap_log3_cmd ln -s $(realpath "${_polap_var_outdir_lk_fq_gz}") "${_polap_var_outdir_nk4_fq_gz}"
	# 	else
	# 		local _RATE=$(echo "scale=4; $limit_file_size / $file_size" | bc)
	# 		# Compare value with 0
	# 		if echo "${_RATE} > 0" | bc -l | grep -q 1; then
	# 			_polap_lib_random-get
	# 			local _random_seed=${_polap_var_random_number}
	# 			_polap_log1 "  random seed for reducing the whole-genome assembly long-read data: ${_random_seed}"
	# 			_polap_log3_pipe "seqkit sample \
	#            -p ${_RATE} \
	#            -s ${_random_seed} \
	#            --quiet \
	#            ${_polap_var_outdir_lk_fq_gz} \
	# 	        -o ${_polap_var_outdir_nk4_fq_gz} \
	#            2>${_polap_output_dest}"
	# 			_polap_log3_pipe "echo ${_random_seed} >${_polap_var_outdir_nk4_fq_gz}.random.seed.${_random_seed}"
	# 		else
	# 			die "ERROR: long-read sampling rate is not greater than 0."
	# 		fi
	# 	fi
	# fi

	_polap_log1 "NEXT (for testing purpose only): $(basename "$0") flye1 -g 150000"
	_polap_log1 "NEXT (for testing purpose only): $(basename "$0") flye1 --test"
	_polap_log1 "NEXT: $(basename "$0") flye1 -o ${_arg_outdir} [-t ${_arg_threads}] [-c ${_arg_coverage}]"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Executes Flye for a whole-genome assembly upto the contigger stage
#
# Arguments:
#   -t ${_arg_threads}: the number of CPU cores
#   -c ${_arg_coverage}: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   ${_arg_outdir}/short_expected_genome_size.txt (ignored with -g option)
#   ${_polap_var_outdir_nk_fq_gz}
# Outputs:
#   $FDIR/30-contigger/contigs.fasta
#   $FDIR/30-contigger/contigs_stats.txt
#   $FDIR/30-contigger/graph_final.fasta
#   $FDIR/30-contigger/graph_final.gfa
################################################################################
function _run_polap_flye1 { # execute Flye for a whole-genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# CHECK: local function
	source "${_POLAPLIB_DIR}/polap-variables-common.sh"
	source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"

	help_message=$(
		cat <<HEREDOC
# Flye whole-genome assembly upto the contigger stage
#
# Arguments:
#   -t ${_arg_threads}: the number of CPU cores
#   --flye-asm-coverage ${_arg_flye_asm_coverage}: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
#   --test
# Inputs:
#   ${_polap_var_outdir_nk_fq_gz}
#   ${_arg_outdir}/short_expected_genome_size.txt (ignored with -g option)
# Outputs:
#   ${_polap_var_wga_contigger_edges_gfa}
#   ${_polap_var_ga_contigger_edges_fasta}
Example: $(basename "$0") ${_arg_menu[0]}
Example: $(basename "$0") ${_arg_menu[0]} --test
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# check_file_existence "${_arg_outdir}/short_expected_genome_size.txt"
	# check_file_existence "${_polap_var_outdir_nk_fq_gz}"

	if [[ -n "${_arg_genomesize}" ]]; then
		_arg_genomesize=$(_polap_utility_convert_unit_to_bp "${_arg_genomesize}")
	fi

	_polap_log0 "flye whole-genome assembly using the reduced long-read data: ${_polap_var_outdir_nk_fq_gz}"
	_polap_log1 "  input1: ${_polap_var_outdir_genome_size}"
	_polap_log1 "  input2: ${_polap_var_outdir_nk_fq_gz}"

	if [ -z "$_arg_genomesize" ]; then
		if [[ ! -s "${_polap_var_outdir_genome_size}" ]]; then
			return ${_POLAP_ERR_NO_GENOME_SIZE}
		fi
	fi

	if [[ "${_arg_long_reads_is}" == "off" ]]; then
		if [[ ! -s "${_polap_var_outdir_nk_fq_gz}" ]]; then
			return ${_POLAP_ERR_NO_NK_FQ}
		fi
	else
		_polap_var_outdir_nk_fq_gz="${_arg_long_reads}"
	fi

	#   ${_polap_var_wga_contigger_edges_gfa}
	if [ -s "${_polap_var_wga_contigger_edges_gfa}" ] &&
		[ "${_arg_redo}" = "off" ]; then
		_polap_log1 "  found: ${_polap_var_wga_contigger_edges_gfa}, so skipping the whole-genome assembly ..."
	else

		if [[ -s "${_polap_var_outdir_genome_size}" ]]; then
			local _EXPECTED_GENOME_SIZE=$(<"${_polap_var_outdir_genome_size}")
			local _EXPECTED_GENOME_SIZE=${_EXPECTED_GENOME_SIZE%.*}
			local _expected_genome_size_bp=$(_polap_utility_convert_bp ${_EXPECTED_GENOME_SIZE})
			_polap_log1 "  expected genome size using short-read data (bases): ${_expected_genome_size_bp}"
		fi
		if [[ ${_arg_test} == "on" ]]; then
			_arg_genomesize=150000
			local _expected_genome_size_bp=$(_polap_utility_convert_bp ${_arg_genomesize})
			_polap_log1 "  test expected genome size (bases): ${_expected_genome_size_bp}"
		fi
		if [ ! -z "$_arg_genomesize" ]; then
			local _EXPECTED_GENOME_SIZE=$_arg_genomesize
			local _expected_genome_size_bp=$(_polap_utility_convert_bp ${_EXPECTED_GENOME_SIZE})
			_polap_log0 "  OPTION: short reads expected genome size (bases) we use instead: ${_expected_genome_size_bp}"
		fi

		_polap_log1 "  executing the whole-genome assembly using flye ... be patient!"
		if [[ "${_arg_timing}" == "off" ]]; then
			local _command1="flye"
		else
			local _command1="command time -v flye"
		fi
		_command1+=" \
      ${_arg_flye_data_type} \
      ${_polap_var_outdir_nk_fq_gz} \
			--out-dir ${_polap_var_wga} \
			--threads ${_arg_threads}"
		if [[ "${_arg_flye_asm_coverage}" -gt 0 ]]; then
			_command1+=" \
			--asm-coverage ${_arg_flye_asm_coverage} \
			--genome-size ${_EXPECTED_GENOME_SIZE}"
		fi
		if [[ "${_arg_menu[1]}" == "polishing" ]]; then
			_command1+=" \
		  --stop-after polishing"
		elif [[ "${_arg_menu[1]}" == "resume" ]]; then
			_command1+=" \
		  --resume"
		else
			_command1+=" \
		  --stop-after contigger"
		fi
		if [[ "${_arg_timing}" == "off" ]]; then
			_command1+=" \
		  2>${_polap_output_dest}"
		else
			_command1+=" \
		  2>${_arg_outdir}/timing-flye1.txt"
		fi
		_polap_log3_pipe "${_command1}"

	fi

	rm -f "${_polap_var_output_wga_gfa}"
	ln -sf "${_polap_var_wga_contigger_edges_gfa}" "${_polap_var_output_wga_gfa}"
	_polap_log1 "  assembly graph in the flye contigger stage: ${_polap_var_wga_contigger_edges_gfa}"
	_polap_log1 "  assembly graph in the flye contigger stage: ${_polap_var_output_wga_gfa}"
	_polap_log1 "NEXT: $(basename "$0") edges-stats -o ${_arg_outdir} [-i ${_arg_inum}]"
	_polap_log1 "NEXT: $(basename "$0") annotate -o ${_arg_outdir} [-i ${_arg_inum}]"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
