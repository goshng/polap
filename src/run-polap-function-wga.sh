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
source "$script_dir/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

source "$script_dir/polap-function-set-variables.sh"
source "$script_dir/run-polap-function-utilities.sh"

################################################################################
# Runs the whole-genome assembly.
################################################################################

function _polap_summary-generated-reads() {
	local _fastq_gz=$1
	local _fastq_stats=$2

	if [ -s "${_fastq_stats}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log1 "  found: ${_fastq_stats}, so skipping POLAP long-read statisics ..."
		_polap_log2_cat "${_fastq_stats}"
	else
		if [ -s "${_fastq_gz}" ]; then
			_polap_log3_pipe "seqkit stats -T ${_fastq_gz} \
				>${_fastq_stats}"
			_polap_log1 "  output2: ${_fastq_stats}"
			_polap_log2_cat "${_fastq_stats}"
		else
			_polap_log1 "  no such file: ${_fastq_gz}, so skipping POLAP long-read statisics ..."
		fi
	fi
}

################################################################################
# Statisics of the short-read and POLAP's long-read (nk.fq.gz) dataset.
################################################################################
function _run_polap_summary-reads() { # statisics of the read dataset
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set paths for bioproject data
	_polap_set-variables-short-read
	source "$script_dir/polap-variables-common.sh"
	source "$script_dir/polap-variables-common.sh"
	source "$script_dir/polap-variables-common.sh"

	help_message=$(
		cat <<HEREDOC
# Summarize the statistics of the long- and short-read dataset.
#
# Arguments:
#   -a ${_arg_short_read1}: a short-read fastq data file
#   -b ${_arg_short_read2}: another short-read fastq data file (optional)
#   or
#   --bioproject use
#   -o ${ODIR}
# Inputs:
#   ${_arg_short_read1}: a short-read fastq data file
#   ${_arg_short_read2}: another short-read fastq data file (optional)
#   ${_polap_var_base_nk_fq_gz}: POLAP generate nk.fq.gz
#   ${_polap_var_base_lk_fq_gz}: POLAP generate lk.fq.gz
# Outputs:
#   ${_polap_var_base_fq_stats}: short-read data statisics
#   ${_polap_var_base_nk_fq_stats}: POLAP long-read data statisics
#   ${_polap_var_base_lk_fq_stats}: POLAP long-read data statisics
# Precondition:
#   (for BioProjectID case)
#   get-bioproject --bioproject <BioProjectID> -o ${ODIR}
Example: $0 ${_arg_menu[0]} -a <file> [-b <file>]
Example: $0 ${_arg_menu[0]} -o ${ODIR} --bioproject use
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [ -s "${_polap_var_base_fq_stats}" ]; then
			_polap_log0_file "${_polap_var_base_fq_stats}"
			_polap_log0_cat "${_polap_var_base_fq_stats}"
		else
			_polap_log0 "No short-read statisics"
		fi

		if [ -s "${_polap_var_base_nk_fq_stats}" ]; then
			_polap_log0_file "${_polap_var_base_nk_fq_stats}"
			_polap_log0_cat "${_polap_var_base_nk_fq_stats}"
		else
			_polap_log0 "No POLAP long-read statisics"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "computing stats of the short-read sequencing data ..."
	if [[ -d "${ODIR}" ]]; then
		_polap_log2 "  output folder: ${ODIR}"
	else
		_polap_log3_cmd mkdir -p "${ODIR}"
	fi
	check_file_existence "${_arg_short_read1}"
	_polap_log1 "  input1: ${_arg_short_read1}"
	if [ -s "${_arg_short_read2}" ]; then
		_polap_log1 "  input2: ${_arg_short_read2}"
	else
		_polap_log1 "  input2: no such file: ${_arg_short_read2}"
		_polap_log1 "    we use a single short-read data file."
	fi

	if [ -s "${_polap_var_base_fq_stats}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log1 "  found: ${_polap_var_base_fq_stats}, so skipping short-read statisics ..."
	else
		if [[ -s "${_arg_short_read1}" ]]; then

			_polap_log3_pipe "seqkit stats -T ${_arg_short_read1} |\
      csvtk del-header >${_polap_var_base_fq_stats}"

			if [ -s "${_arg_short_read2}" ]; then
				_polap_log3_pipe "seqkit stats -T ${_arg_short_read2} |\
				csvtk del-header >>${_polap_var_base_fq_stats}"
			fi
		else
			_polap_log0 "  no such file: ${_arg_short_read1}"
		fi
	fi

	_polap_log1 "  output1: ${_polap_var_base_fq_stats}"
	_polap_log2_cat "${_polap_var_base_fq_stats}"

	_polap_summary-generated-reads \
		"${_polap_var_base_nk_fq_gz}" \
		"${_polap_var_base_nk_fq_stats}"

	_polap_summary-generated-reads \
		"${_polap_var_base_lk_fq_gz}" \
		"${_polap_var_base_lk_fq_stats}"

	_polap_log1 NEXT: $0 total-length-long -o "$ODIR" -l ${_arg_long_reads}
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
#   $ODIR
# Outputs:
#   $ODIR/long_total_length.txt
################################################################################
function _run_polap_total-length-long() { # total size (bp) of long-read data
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# CHECK: local function
	_polap_set-variables-long-read
	source "$script_dir/polap-variables-common.sh"
	source "$script_dir/run-polap-function-utilities.sh"

	help_message=$(
		cat <<HEREDOC
# Compute the total number of nucleotides of long-read data.
#
# Arguments:
#   -l ${_arg_long_reads}: a long-read fastq data file (the highest priority)
#   or
#   --bioproject use
#   -o ${ODIR}: ${ODIR}/0-bioproject (the least priority)
# Inputs:
#   ${_arg_long_reads}: a long-read fastq data file
# Outputs:
#   ${_polap_var_base_long_total_length}
# Precondition:
#   (for BioProjectID case)
#   get-bioproject --bioproject <BioProjectID> -o ${ODIR}
Example: $(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads}
Example: $(basename "$0") ${_arg_menu[0]} -o ${ODIR} --bioproject use
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [ -s "${_polap_var_base_long_total_length}" ]; then
			_polap_log0_cat "${_polap_var_base_long_total_length}"
		else
			_polap_log0 "No long-read total length"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		return
	fi

	_polap_log0 "counting the total number of bases in the long-read dataset of $_arg_long_reads ..."
	if [[ ! -d "${ODIR}" ]]; then
		_polap_log2 "  creating output folder: ${ODIR}"
		_polap_log3_cmd mkdir -p "${ODIR}"
	fi
	check_file_existence "${_arg_long_reads}"

	_polap_log1 "  input: ${_arg_long_reads}"

	if [ -s "${_polap_var_base_long_total_length}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log1 "  found: ${_polap_var_base_long_total_length}, so skipping long-read statisics ..."
	else
		_polap_log3_pipe "seqkit stats -Ta ${_arg_long_reads} |\
			csvtk cut -t -f sum_len |\
			csvtk del-header \
				>${_polap_var_base_long_total_length}"
	fi

	_polap_log1 "  output: ${_polap_var_base_long_total_length}"
	_polap_log2_cat "${_polap_var_base_long_total_length}"

	local LONG_TOTAL_LENGTH=$(<"${_polap_var_base_long_total_length}")
	local _total_long_read=$(_polap_utility_convert_bp ${LONG_TOTAL_LENGTH})
	_polap_log0 "  total length of the long-read dataset (bases): ${_total_long_read}"

	_polap_log1 NEXT: $0 find-genome-size -o "$ODIR" -a "${_arg_short_read1}" -b "${_arg_short_read2}"
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
#   $ODIR/jellyfish_out
#   $ODIR/jellyfish_out.histo
#   $ODIR/short_expected_genome_size.txt
################################################################################
function _run_polap_find-genome-size() { # estimate the whole genome size
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Set paths for bioproject data
	_polap_set-variables-short-read
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'
	source "$script_dir/polap-variables-common.sh" # '.' means 'source'
	source "$script_dir/run-polap-function-utilities.sh"

	help_message=$(
		cat <<HEREDOC
# Estimates the whole genome size using short-read data.
#
# Arguments:
#   -a ${_arg_short_read1}: a short-read fastq data file
#   -b ${_arg_short_read2}: another short-read fastq data file (optional)
#   or
#   --bioproject use
#   -o ${ODIR}
# Inputs:
#   ${_arg_short_read1}: a short-read fastq data file
#   ${_arg_short_read2}: another short-read fastq data file (optional)
# Outputs:
#   ${_polap_var_base_genome_size}
Example: $(basename "$0") ${_arg_menu[0]} -a <file> [-b <file>]
Example: $(basename "$0") ${_arg_menu[0]} -o ${ODIR} --bioproject use
HEREDOC
	)

	# Display help message
	[[ "${_arg_menu[1]}" == "help" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [ -s "${_polap_var_base_genome_size}" ]; then
			_polap_log0_file "${_polap_var_base_genome_size}"
			_polap_log0_cat "${_polap_var_base_genome_size}"
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

	if [[ -d "${ODIR}" ]]; then
		_polap_log2 "  output folder: ${ODIR}"
	else
		_polap_log3 mkdir -p "${ODIR}"
		mkdir -p "${ODIR}"
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
	if [ -s "${_polap_var_base_genome_size}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log1 "  found: ${_polap_var_base_genome_size}, so skipping the genome size estimation using the short-read data ..."
		_polap_log1_file "${_polap_var_base_genome_size}"
	else
		if [ -s "${_polap_var_base_jellyfish_out}" ] && [ "${_arg_redo}" = "off" ]; then
			_polap_log1 "  found: ${_polap_var_base_jellyfish_out}, so skipping the JellyFish count step using the short-read data ..."
			_polap_log1_file "${_polap_var_base_jellyfish_out}"
		else
			if [[ -s "${_arg_short_read1}" ]] && [[ -s "${_arg_short_read2}" ]]; then
				_polap_log3_pipe "jellyfish count \
					-t ${_arg_threads} -C -m 19 \
					-s 5G \
					-o ${_polap_var_base_jellyfish_out} \
					--min-qual-char=? \
					${_arg_short_read1} ${_arg_short_read2}"
			elif [[ -s "${_arg_short_read1}" ]] && [[ ! -s "${_arg_short_read2}" ]]; then
				_polap_log3_pipe "jellyfish count \
					-t ${_arg_threads} -C -m 19 \
					-s 5G \
					-o ${_polap_var_base_jellyfish_out} \
					--min-qual-char=? \
					${_arg_short_read1}"
			else
				die "ASSERT: we must have at least one short-read fastq file."
			fi
			check_file_existence "${_polap_var_base_jellyfish_out}"
		fi

		check_file_existence "${_polap_var_base_jellyfish_out}"
		_polap_log3_pipe "jellyfish histo \
			-o ${_polap_var_base_jellyfish_out_histo} \
			${_polap_var_base_jellyfish_out}"
		_polap_log2_file "${_polap_var_base_jellyfish_out_histo}"

		_polap_log3_pipe "Rscript $script_dir/run-polap-r-jellyfish.R \
			${_polap_var_base_jellyfish_out_histo} \
			${_polap_var_base_genome_size}"
		# Check the exit status
		if [ $? -ne 0 ]; then
			# Take action if needed, e.g., logging, sending a notification, etc.
			die "ERROR: not being to estimate the genome size using short-read data: short-read may be too small."
		fi
	fi
	_polap_log1 "  output: ${_polap_var_base_genome_size}"
	_polap_log2_cat "${_polap_var_base_genome_size}"

	EXPECTED_GENOME_SIZE=$(<"${_polap_var_base_genome_size}")
	EXPECTED_GENOME_SIZE=${EXPECTED_GENOME_SIZE%.*}
	local _expected_genome_size=$(_polap_utility_convert_bp ${EXPECTED_GENOME_SIZE})
	_polap_log0 "  expected genome size using short-read data (bases): ${_expected_genome_size}"

	_polap_log1 NEXT: $0 reduce-data -o "$ODIR" -l "${_arg_long_reads}" [-m "${_arg_min_read_length}"]

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

################################################################################
# Checks if the long-read coverage is less than ${_arg_coverage}.
#
# If so, keep the long read data.
# If not, sample long reads upto that coverage.
# Deletes long reads shorter than a sequence length threshold e.g., 3 kb.
#
# Arguments:
#   -l $_arg_long_reads: a long-read fastq data file
#   -m ${_arg_min_read_length}: the long-read sequence length threshold
#   --reduction-reads (default) or --no-reduction-reads
# Inputs:
#   $ODIR/short_expected_genome_size.txt
#   $ODIR/long_total_length.txt
#   $_arg_long_reads
# Outputs:
#   ${_polap_var_base_nk_fq_gz}
################################################################################
function _run_polap_reduce-data() { # reduce the long-read data, if too big
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# CHECK: local function
	source "$script_dir/polap-variables-common.sh"
	_polap_set-variables-long-read
	source "$script_dir/run-polap-function-utilities.sh"

	help_message=$(
		cat <<HEREDOC
# FIXME: test it. Reduce the long-read data for a Flye genome assembly.
# Reduce the long-read data for a Flye genome assembly.
#
# Outputs:
# 1. ${_polap_var_base_lk_fq_gz}: reads longer than ${_arg_min_read_length} bp
# 2. ${_polap_var_base_nk_fq_gz}: subsample of ${_arg_long_reads}
#
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
#
# Arguments:
#   -l $_arg_long_reads: a long-read fastq data file
#   or
#   --bioproject use
#   -o ${ODIR}: ${ODIR}/0-bioproject (the least priority)
#   -m ${_arg_min_read_length}: the long-read sequence length threshold
#   -c ${_arg_coverage}: the target coverage
#   --no-reduction-reads
#   --random-seed 11: for seqkit default seed
#   --test
# Inputs:
#   ${_polap_var_base_genome_size}
#   ${_polap_var_base_long_total_length}
#   ${_arg_long_reads}
# Outputs:
#   ${_polap_var_base_lk_fq_gz}: reads longer than ${_arg_min_read_length} bp
#   ${_polap_var_base_nk_fq_gz}: subsample of ${_polap_var_base_nk_fq_gz}
# Menu:
#   split <N>: splits ${_polap_var_base_lk_fq_gz} into N parts.
Example: $0 ${_arg_menu[0]} -l <arg> -m <arg>
Example: $0 ${_arg_menu[0]} -o ${ODIR} --bioproject use
Example: $0 ${_arg_menu[0]} split 10
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [ -s "${_polap_var_base_nk_fq_gz}" ]; then
			_polap_log0_cat "${_polap_var_base_nk_fq_stats}"
		else
			_polap_log0 "No reduced subsample long-read data file: ${_polap_var_base_nk_fq_gz}"
		fi

		if [ -s "${_polap_var_base_lk_fq_gz}" ]; then
			_polap_log0_cat "${_polap_var_base_lk_fq_stats}"
		else
			_polap_log0 "No reduced long-read: ${_polap_var_base_lk_fq_gz}"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	if [[ "${_arg_menu[1]}" == "split" ]]; then
		if [ -s "${_polap_var_base_lk_fq_gz}" ]; then
			seqkit split2 -p "${_arg_menu[2]}" "${_polap_var_base_lk_fq_gz}"
		else
			_polap_log0 "No reduced long-read: ${_polap_var_base_lk_fq_gz}"
		fi

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "reducing the long-read data ${_arg_long_reads} for the whole- and organelle-genome assemblies ..."

	# Check for required files
	check_folder_existence "${ODIR}"
	check_file_existence "${_polap_var_base_long_total_length}"
	check_file_existence "${_polap_var_base_genome_size}"
	check_file_existence "${_arg_long_reads}"

	_polap_log1 "  input1: ${_arg_long_reads}"
	_polap_log1 "  input2: ${_polap_var_base_long_total_length}"
	_polap_log1 "  input3: ${_polap_var_base_genome_size}"

	# Get the expected genome size and long-read sequencing coverage
	local _LONG_TOTAL_LENGTH=$(<"${_polap_var_base_long_total_length}")
	local _long_total_length_bp=$(_polap_utility_convert_bp ${_LONG_TOTAL_LENGTH})
	local _EXPECTED_GENOME_SIZE=$(<"${_polap_var_base_genome_size}")
	local _expected_genome_size_bp=$(_polap_utility_convert_bp ${_EXPECTED_GENOME_SIZE})
	local _EXPECTED_LONG_COVERAGE=$(echo "scale=3; ${_LONG_TOTAL_LENGTH}/${_EXPECTED_GENOME_SIZE}" | bc)
	_polap_log2 "    LONG_TOTAL_LENGTH: ${_long_total_length_bp}"
	_polap_log2 "    EXPECTED_GENOME_SIZE: ${_expected_genome_size_bp}"
	_polap_log2 "    EXPECTED_LONG_COVERAGE: ${_EXPECTED_LONG_COVERAGE}x"

	if [[ -s "${_polap_var_base_nk_fq_gz}" ]] &&
		[[ -s "${_polap_var_base_lk_fq_gz}" ]] &&
		[[ "${_arg_redo}" == "off" ]]; then
		_polap_log0 "  found: ${_polap_var_base_nk_fq_gz}, so skipping the long-read data reduction."
		_polap_log0 "  found: ${_polap_var_base_lk_fq_gz}, so skipping the long-read data reduction."
		return
	fi

	_polap_log0 "  reducing ${_arg_long_reads} -> ${_polap_var_base_lk_fq_gz} with the minimum length of ${_arg_min_read_length} bp"

	if [[ -s "${_polap_var_base_lk_fq_gz}" ]] && [[ "${_arg_redo}" == "off" ]]; then
		_polap_log0 "  found: ${_polap_var_base_lk_fq_gz}, so skipping the long-read data reduction."
	else
		_polap_log1 "  keeps long reads of length being at least ${_arg_min_read_length} bp ..."
		_polap_log2 "    input1: ${_arg_long_reads}"
		_polap_log2 "    output: ${_polap_var_base_lk_fq_gz}"
		_polap_log2 "    sorting reads for a faster seqtk subseq ..."
		_polap_log3_cmd rm -f "${_polap_var_base_lk_fq_gz}"
		_polap_log3_pipe "seqkit seq \
      --quiet \
      -m ${_arg_min_read_length} \
      --threads 4 \
		  ${_arg_long_reads} \
		  -o ${_polap_var_base_lk_fq_gz} \
      >${_polap_output_dest} 2>&1"

		_polap_log1 "  creating the statisics for the size reduced long-read data ..."
		_polap_log2 "    input1: ${_polap_var_base_lk_fq_gz}"
		_polap_log2 "    output: ${_polap_var_base_lk_fq_stats}"
		_polap_log3_pipe "seqkit stats -aT \
      ${_polap_var_base_lk_fq_gz} \
      >${_polap_var_base_lk_fq_stats}"
		_polap_log2_column "${_polap_var_base_lk_fq_stats}"
	fi

	_polap_log0 "  reducing ${_arg_long_reads} -> ${_polap_var_base_nk_fq_gz} with the target coverage of ${_arg_coverage}x ..."

	if [[ -s "${_polap_var_base_nk_fq_gz}" ]] && [[ "${_arg_redo}" == "off" ]]; then
		_polap_log0 "  found: ${_polap_var_base_nk_fq_gz}, so skipping the long-read data reduction."
	else
		# subsample the long-read data so that the target coverage is ${_arg_coverage}.
		local nfq_file="${ODIR}/n.fq"
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
				local _RATE=$(echo "scale=3; ${_arg_coverage}/${_EXPECTED_LONG_COVERAGE}" | bc)
				# Compare value with 0
				if echo "${_RATE} > 0" | bc -l | grep -q 1; then
					_polap_log1 "  sampling long-read data by $_RATE ..."
					_polap_log1 "    $_RATE <= target long-read genome coverage[${_arg_coverage}]/expected long-read genome coverage[$EXPECTED_LONG_COVERAGE] ..."
					local _random_seed=${_arg_random_seed:-$RANDOM}
					_polap_log0 "  random seed for reducing the whole-genome assembly long-read data: ${_random_seed}"
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
		          -o ${_polap_var_base_nk_fq_gz} \
              >${_polap_output_dest} 2>&1"
					_polap_log3_cmd touch "${_polap_var_base_nk_fq_gz}.random.seed.${_random_seed}"
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
			_polap_log2 "    output: ${_polap_var_base_nk_fq_gz}"
			_polap_log3_cmd rm -f "${_polap_var_base_nk_fq_gz}"
			_polap_log3_pipe "seqkit seq \
        --quiet \
        -m ${_arg_min_read_length} \
        --threads 4 \
		    ${nfq_file} \
  		  -o ${_polap_var_base_nk_fq_gz} \
        >${_polap_output_dest} 2>&1"

			_polap_log3_cmd rm -f "$nfq_file"
		fi

		_polap_log1 "  creating the statisics for the reduced long-read data ..."
		_polap_log3_pipe "seqkit stats -aT \
      ${_polap_var_base_nk_fq_gz} \
      >${_polap_var_base_nk_fq_stats}"
		_polap_log1 "    output: ${_polap_var_base_nk_fq_stats}"
		_polap_log2_column "${_polap_var_base_nk_fq_stats}"
	fi

	_polap_log1 "NEXT (for testing purpose only): $0 flye1 -g 150000"
	_polap_log1 "NEXT (for testing purpose only): $0 flye1 --test"
	_polap_log1 "NEXT: $0 flye1 -o $ODIR [-t ${_arg_threads}] [-c ${_arg_coverage}]"

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
#   $ODIR/short_expected_genome_size.txt (ignored with -g option)
#   ${_polap_var_base_nk_fq_gz}
# Outputs:
#   $FDIR/30-contigger/contigs.fasta
#   $FDIR/30-contigger/contigs_stats.txt
#   $FDIR/30-contigger/graph_final.fasta
#   $FDIR/30-contigger/graph_final.gfa
################################################################################
function _run_polap_flye1() { # execute Flye for a whole-genome assembly
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# CHECK: local function
	source "$script_dir/polap-variables-common.sh"
	source "$script_dir/polap-variables-common.sh"
	source "$script_dir/run-polap-function-utilities.sh"

	LRNK="${_polap_var_base_nk_fq_gz}"

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
#   ${_polap_var_base_nk_fq_gz}
#   $ODIR/short_expected_genome_size.txt (ignored with -g option)
# Outputs:
#   ${_polap_var_wga_contigger_gfa}
#   ${_polap_var_contigger_edges_fasta}
#   ${_polap_var_wga_contigger_contigs_stats}
#   ${_polap_var_wga_contigger_contigs_fasta}
Example: $0 ${_arg_menu[0]}
Example: $0 ${_arg_menu[0]} --test
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	check_file_existence "$ODIR/short_expected_genome_size.txt"
	check_file_existence "${_polap_var_base_nk_fq_gz}"

	_polap_log0 "flye whole-genome assembly using the reduced long-read data: ${_polap_var_base_nk_fq_gz}"
	_polap_log1 "  input1: ${_polap_var_base_genome_size}"
	_polap_log1 "  input2: ${_polap_var_base_nk_fq_gz}"

	#   ${_polap_var_wga_contigger_gfa}
	if [ -s "${_polap_var_wga_contigger_gfa}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log1 "  found: ${_polap_var_wga_contigger_gfa}, so skipping the whole-genome assembly ..."
	else
		local EXPECTED_GENOME_SIZE=$(<"${_polap_var_base_genome_size}")
		local EXPECTED_GENOME_SIZE=${EXPECTED_GENOME_SIZE%.*}
		local _expected_genome_size=$(_polap_utility_convert_bp ${EXPECTED_GENOME_SIZE})
		_polap_log1 "  expected genome size using short-read data (bases): ${_expected_genome_size}"
		if [[ ${_arg_test} == "on" ]]; then
			_arg_genomesize=150000
		fi
		if [ ! -z "$_arg_genomesize" ]; then
			EXPECTED_GENOME_SIZE=$_arg_genomesize
			local _expected_genome_size=$(_polap_utility_convert_bp ${EXPECTED_GENOME_SIZE})
			_polap_log0 "  OPTION: short reads expected genome size (bases) we use instead: ${_expected_genome_size}"
		fi

		_polap_log1 "  executing the whole-genome assembly using flye ... be patient!"
		# _polap_log3_pipe "flye \
		#     --nano-raw ${_polap_var_base_nk_fq_gz} \
		# 	--out-dir ${_polap_var_wga} \
		# 	--threads ${_arg_threads} \
		# 	--asm-coverage ${_arg_flye_asm_coverage} \
		# 	--genome-size ${EXPECTED_GENOME_SIZE} \
		# 	--stop-after contigger \
		# 	>${_polap_output_dest} 2>&1"

		local _command1="flye \
      --nano-raw ${_polap_var_base_nk_fq_gz} \
			--out-dir ${_polap_var_wga} \
			--threads ${_arg_threads}"
		if [[ "${_arg_flye_asm_coverage}" -gt 0 ]]; then
			_command1+=" \
			--asm-coverage ${_arg_flye_asm_coverage} \
			--genome-size ${EXPECTED_GENOME_SIZE}"
		fi
		if [[ "${_arg_menu[2]}" == "polishing" ]]; then
			_command1+=" \
		  --resume"
		else
			_command1+=" \
		  --stop-after contigger"
		fi
		_command1+=" \
		  2>${_polap_output_dest}"
		_polap_log3_pipe "${_command1}"

	fi

	_polap_log1 "  assembly graph in the flye contigger stage: ${_polap_var_wga_contigger_gfa}"
	_polap_log1 "NEXT: $0 edges-stats -o $ODIR [-i ${INUM}]"
	_polap_log1 "NEXT: $0 annotate -o $ODIR [-i ${INUM}]"

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}
