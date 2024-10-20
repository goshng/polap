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
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

source "$script_dir/polap-function-set-variables.sh"

################################################################################
# Runs the whole-genome assembly.
################################################################################

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
	source "$script_dir/polap-variables-base.sh"
	source "$script_dir/polap-variables-bioproject.sh"
	source "$script_dir/polap-variables-oga.sh"

	help_message=$(
		cat <<HEREDOC
# Summarize the statistics of the long- and short-read dataset.
#
# Arguments:
#   -a $SR1: a short-read fastq data file
#   -b $SR2: another short-read fastq data file (optional)
#   or
#   --bioproject use
#   -o ${ODIR}
# Inputs:
#   $SR1: a short-read fastq data file
#   $SR2: another short-read fastq data file (optional)
# Outputs:
#   ${_polap_var_base_fq_stats}: short-read data statisics
#   ${_polap_var_base_nk_fq_stats}: POLAP long-read data statisics
# Precondition:
#   (for BioProjectID case)
#   get-bioproject --bioproject <BioProjectID> -o ${ODIR}
Example: $(basename "$0") ${_arg_menu[0]} -a <file> [-b <file>]
Example: $(basename "$0") ${_arg_menu[0]} -o ${ODIR} --bioproject use
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

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

		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "computing stats of the short-read sequencing data ..."
	if [[ -d "${ODIR}" ]]; then
		_polap_log2 "  output folder: ${ODIR}"
	else
		_polap_log3 mkdir -p "${ODIR}"
		mkdir -p "${ODIR}"
	fi
	check_file_existence "${SR1}"
	_polap_log0 "  input1: ${SR1}"
	if [ -s "${SR2}" ]; then
		_polap_log0 "  input2: ${SR2}"
	else
		_polap_log0 "  input2: no such file: ${SR2}"
		_polap_log0 "    we use a single short-read data file."
	fi

	if [ -s "${_polap_var_base_fq_stats}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log0 "  found: ${_polap_var_base_fq_stats}, so skipping short-read statisics ..."
	else
		_polap_log3 seqkit stats -T "$SR1" "|" csvtk del-header ">${_polap_var_base_fq_stats}"
		seqkit stats -T "$SR1" |
			csvtk del-header >"${_polap_var_base_fq_stats}"

		if [ -s "${SR2}" ]; then
			_polap_log3 seqkit stats -T "$SR2" "|" csvtk del-header ">>${_polap_var_base_fq_stats}"
			seqkit stats -T "$SR2" |
				csvtk del-header >>"${_polap_var_base_fq_stats}"
		fi
	fi

	_polap_log0 "  output1: ${_polap_var_base_fq_stats}"
	_polap_log2_cat "${_polap_var_base_fq_stats}"

	if [ -s "${_polap_var_base_nk_fq_stats}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log0 "  found: ${_polap_var_base_nk_fq_stats}, so skipping POLAP long-read statisics ..."
		_polap_log2_cat "${_polap_var_base_nk_fq_stats}"
	else
		if [ -s "${_polap_var_base_nk_fq_gz}" ]; then
			_polap_log3 seqkit stats -T "${_polap_var_base_nk_fq_gz}" ">${_polap_var_base_nk_fq_stats}"
			seqkit stats -T "${_polap_var_base_nk_fq_gz}" \
				>"${_polap_var_base_nk_fq_stats}"
			_polap_log0 "${_polap_var_base_nk_fq_stats}"
			_polap_log2_cat "${_polap_var_base_nk_fq_stats}"
		else
			_polap_log0 "  no such file: ${_polap_var_base_nk_fq_gz}, so skipping POLAP long-read statisics ..."
		fi
	fi

	_polap_log1 NEXT: $(basename "$0") total-length-long -o "$ODIR" -l ${_arg_long_reads}
	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Computes the total number of nucleotides of long-read data.
#
# Arguments:
#   -l $LR: a long-read fastq data file
# Inputs:
#   $LR: a long-read fastq data file
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
	source "$script_dir/polap-variables-base.sh"
	source "$script_dir/polap-variables-bioproject.sh"
	source "$script_dir/polap-variables-oga.sh"
	source "$script_dir/run-polap-function-utilities.sh"

	help_message=$(
		cat <<HEREDOC
# Computes the total number of nucleotides of long-read data.
#
# Arguments:
#   -l ${LR}: a long-read fastq data file (the highest priority)
#   or
#   --bioproject use
#   -o ${ODIR}: ${ODIR}/0-bioproject (the least priority)
# Inputs:
#   ${LR}: a long-read fastq data file
# Outputs:
#   ${_polap_var_base_long_total_length}
# Precondition:
#   (for BioProjectID case)
#   get-bioproject --bioproject <BioProjectID> -o ${ODIR}
Example: $(basename "$0") ${_arg_menu[0]} [-l|--long-reads ${LR}]
Example: $(basename "$0") ${_arg_menu[0]} -o ${ODIR} --bioproject use
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [ -s "${_polap_var_base_long_total_length}" ]; then
			_polap_log0_file "${_polap_var_base_long_total_length}"
			_polap_log0_cat "${_polap_var_base_long_total_length}"
		else
			_polap_log0 "No long-read total length"
		fi

		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "counting the total number of bases in the long-read dataset of $LR ..."
	if [[ -d "${ODIR}" ]]; then
		_polap_log2 "  output folder: ${ODIR}"
	else
		_polap_log3 mkdir -p "${ODIR}"
		mkdir -p "${ODIR}"
	fi
	check_file_existence "${LR}"

	_polap_log0 "  input: ${LR}"

	if [ -s "${_polap_var_base_long_total_length}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log0 "  found: ${_polap_var_base_long_total_length}, so skipping long-read statisics ..."
	else
		_polap_log3 seqkit stats -Ta "${LR}" "|" csvtk cut -t -f "sum_len" "|" csvtk del-header ">${_polap_var_base_long_total_length}"
		seqkit stats -Ta "${LR}" |
			csvtk cut -t -f "sum_len" |
			csvtk del-header \
				>"${_polap_var_base_long_total_length}"
	fi

	_polap_log0 "  output: ${_polap_var_base_long_total_length}"
	_polap_log2_cat "${_polap_var_base_long_total_length}"

	local LONG_TOTAL_LENGTH=$(<"${_polap_var_base_long_total_length}")
	local _total_long_read=$(_polap_utility_convert_bp ${LONG_TOTAL_LENGTH})
	_polap_log0 "  total length of the long-read dataset (bases): ${_total_long_read}"

	_polap_log1 NEXT: $(basename "$0") find-genome-size -o "$ODIR" -a "${_arg_short_read1}" -b "${_arg_short_read2}"
	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
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
#   -a $SR1: a short-read fastq data file
#   -b $SR2: another short-read fastq data file (optional)
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
	source "$script_dir/polap-variables-base.sh"       # '.' means 'source'
	source "$script_dir/polap-variables-bioproject.sh" # '.' means 'source'
	source "$script_dir/polap-variables-oga.sh"        # '.' means 'source'
	source "$script_dir/run-polap-function-utilities.sh"

	help_message=$(
		cat <<HEREDOC
# Estimates the whole genome size using short-read data.
#
# Arguments:
#   -a $SR1: a short-read fastq data file
#   -b $SR2: another short-read fastq data file (optional)
#   or
#   --bioproject use
#   -o ${ODIR}
# Inputs:
#   $SR1: a short-read fastq data file
#   $SR2: another short-read fastq data file (optional)
# Outputs:
#   ${_polap_var_base_genome_size}
Example: $(basename "$0") ${_arg_menu[0]} -a <file> [-b <file>]
Example: $(basename "$0") ${_arg_menu[0]} -o ${ODIR} --bioproject use
HEREDOC
	)

	# Display help message
	[[ "${_arg_menu[1]}" == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [ -s "${_polap_var_base_genome_size}" ]; then
			_polap_log0_file "${_polap_var_base_genome_size}"
			_polap_log0_cat "${_polap_var_base_genome_size}"
		else
			_polap_log0 "No genome size estimate from a short-read dataset"
		fi

		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "estimating the genome size using your short-read data [$SR1] and [$SR2] ..."

	if [[ -d "${ODIR}" ]]; then
		_polap_log2 "  output folder: ${ODIR}"
	else
		_polap_log3 mkdir -p "${ODIR}"
		mkdir -p "${ODIR}"
	fi
	check_file_existence "${SR1}"

	_polap_log0 "  input1: $SR1"
	if [ -s "${SR2}" ]; then
		_polap_log0 "  input2: ${SR2}"
	else
		_polap_log0 "  input2: no such file: ${SR2}"
		_polap_log0 "    we use a single short-read data file."
	fi

	# See https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/
	if [ -s "${_polap_var_base_genome_size}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log0 "  found: ${_polap_var_base_genome_size}, so skipping the genome size estimation using the short-read data ..."
		_polap_log1_file "${_polap_var_base_genome_size}"
	else
		if [ -s "${_polap_var_base_jellyfish_out}" ] && [ "${_arg_redo}" = "off" ]; then
			_polap_log0 "  found: ${_polap_var_base_jellyfish_out}, so skipping the JellyFish count step using the short-read data ..."
			_polap_log1_file "${_polap_var_base_jellyfish_out}"
		else
			if [[ -s "$SR1" ]] && [[ -s "$SR2" ]]; then
				_polap_log3 jellyfish count -t "${NT}" -C -m 19 -s 5G -o "${_polap_var_base_jellyfish_out}" --min-qual-char=? "$SR1" "$SR2"
				jellyfish count \
					-t "$NT" -C -m 19 \
					-s 5G \
					-o "${_polap_var_base_jellyfish_out}" \
					--min-qual-char=? \
					"$SR1" "$SR2"
			elif [[ -s "$SR1" ]] && [[ ! -s "$SR2" ]]; then
				_polap_log3 jellyfish count -t "$NT" -C -m 19 -s 5G -o "${_polap_var_base_jellyfish_out}" --min-qual-char=? "$SR1"
				jellyfish count \
					-t "$NT" -C -m 19 \
					-s 5G \
					-o "${_polap_var_base_jellyfish_out}" \
					--min-qual-char=? \
					"$SR1"
			else
				die "ASSERT: we must have at least one short-read fastq file."
			fi
			check_file_existence "${_polap_var_base_jellyfish_out}"
		fi

		check_file_existence "${_polap_var_base_jellyfish_out}"
		_polap_log3 jellyfish histo -o "${_polap_var_base_jellyfish_out_histo}" "${_polap_var_base_jellyfish_out}"
		jellyfish histo \
			-o "${_polap_var_base_jellyfish_out_histo}" \
			"${_polap_var_base_jellyfish_out}"
		_polap_log2_file "${_polap_var_base_jellyfish_out_histo}"

		_polap_log3 Rscript "$script_dir/run-polap-jellyfish.R" \
			"${_polap_var_base_jellyfish_out_histo}" \
			"${_polap_var_base_genome_size}"
		Rscript "$script_dir/run-polap-jellyfish.R" \
			"${_polap_var_base_jellyfish_out_histo}" \
			"${_polap_var_base_genome_size}"
		# Check the exit status
		if [ $? -ne 0 ]; then
			# Take action if needed, e.g., logging, sending a notification, etc.
			die "ERROR: not being to estimate the genome size using short-read data: short-read may be too small."
		fi
	fi
	_polap_log0 "  output: ${_polap_var_base_genome_size}"
	_polap_log1_cat "${_polap_var_base_genome_size}"

	EXPECTED_GENOME_SIZE=$(<"${_polap_var_base_genome_size}")
	EXPECTED_GENOME_SIZE=${EXPECTED_GENOME_SIZE%.*}
	local _expected_genome_size=$(_polap_utility_convert_bp ${EXPECTED_GENOME_SIZE})
	_polap_log0 "  expected genome size using short-read data (bases): ${_expected_genome_size}"

	_polap_log1 NEXT: $(basename "$0") reduce-data -o "$ODIR" -l "${_arg_long_reads}" [-m "${_arg_min_read_length}"]

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Checks if the long-read coverage is less than $COV.
#
# If so, keep the long read data.
# If not, sample long reads upto that coverage.
# Deletes long reads shorter than a sequence length threshold e.g., 3 kb.
#
# Arguments:
#   -l $LR: a long-read fastq data file
#   -m $MR: the long-read sequence length threshold
#   --reduction-reads (default) or --no-reduction-reads
# Inputs:
#   $ODIR/short_expected_genome_size.txt
#   $ODIR/long_total_length.txt
#   $LR
# Outputs:
#   $LRNK
################################################################################
function _run_polap_reduce-data() { # reduce the long-read data, if too big
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# CHECK: local function
	_polap_set-variables-long-read
	source "$script_dir/polap-variables-base.sh"
	source "$script_dir/polap-variables-bioproject.sh"
	source "$script_dir/polap-variables-oga.sh"
	source "$script_dir/run-polap-function-utilities.sh"

	LRNK="${_polap_var_base_nk_fq_gz}"

	help_message=$(
		cat <<HEREDOC
# Reduce the long-read data for the flye genome assembly.
#
# 1. Subsample the long-read data size with a target coverage.
# Checks if the long-read coverage is less than $COV.
# If so, keep the long read data. Sample long reads upto that coverage, otherwise.
# 2. Deletes long reads shorter than a sequence length threshold e.g., 3 kb.
#
# Arguments:
#   -l $LR: a long-read fastq data file
#   or
#   --bioproject use
#   -o ${ODIR}: ${ODIR}/0-bioproject (the least priority)
#
#   -m $MR: the long-read sequence length threshold
#   -c $COV: the target coverage
#   --reduction-reads (default) or --no-reduction-reads
# Inputs:
#   ${_polap_var_base_long_total_length}
#   ${_polap_var_base_genome_size}
#   ${LR}
# Outputs:
#   ${_polap_var_base_nk_fq_gz}
Example: $(basename "$0") ${_arg_menu[0]} -l <arg> -m <arg>
Example: $(basename "$0") ${_arg_menu[0]} -o ${ODIR} --bioproject use
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		if [ -s "${_polap_var_base_nk_fq_gz}" ]; then
			_polap_log0_file "${_polap_var_base_nk_fq_gz}"
			_polap_log0_cat "${_polap_var_base_nk_fq_stats}"
		else
			_polap_log0 "No reduced long-read"
		fi

		_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$DEBUG" -eq 1 ] && set +x
		exit $EXIT_SUCCESS
	fi

	_polap_log0 "reducing the long-read data ${LR} with the target coverage of ${COV}x ..."

	# Check for required files
	check_folder_existence "${ODIR}"
	check_file_existence "${_polap_var_base_long_total_length}"
	check_file_existence "${_polap_var_base_genome_size}"
	check_file_existence "${LR}"

	_polap_log0 "  input1: ${_polap_var_base_long_total_length}"
	_polap_log0 "  input2: ${_polap_var_base_genome_size}"

	# Get the expected genome size and long-read sequencing coverage
	local LONG_TOTAL_LENGTH=$(<"${_polap_var_base_long_total_length}")
	local EXPECTED_GENOME_SIZE=$(<"${_polap_var_base_genome_size}")
	local EXPECTED_GENOME_SIZE=${EXPECTED_GENOME_SIZE%.*}
	local EXPECTED_LONG_COVERAGE=$(echo "scale=3; $LONG_TOTAL_LENGTH/$EXPECTED_GENOME_SIZE" | bc)

	if [[ -s "${LRNK}" ]] && [ "${_arg_redo}" = "off" ]; then
		_polap_log0 "  found: ${LRNK}, so skipping the long-read data reduction."
		[ "$DEBUG" -eq 1 ] && set +x
		return $EXIT_SUCCESS
	fi

	# subsample the long-read data so that the target coverage is $COV.
	local nfq_file="${ODIR}/n.fq"
	_polap_log2 "  deletes ${nfq_file} if there is one."
	_polap_log3 rm -f "${nfq_file}"
	rm -f "${nfq_file}"
	if [[ ${_arg_test} == "on" ]]; then
		_polap_log1 "  OPTION: --test : No reduction of the test long-read data"
		_polap_log3 ln -s $(realpath "$LR") "$nfq_file"
		ln -s $(realpath "$LR") "$nfq_file"
	elif [[ ${_arg_reduction_reads} == "off" ]]; then
		_polap_log1 "  OPTION: --no-reduction-reads : No reduction of the long-read data"
		_polap_log3 ln -s $(realpath "$LR") "$nfq_file"
		ln -s $(realpath "$LR") "$nfq_file"
	else
		# if [ "$EXPECTED_LONG_COVERAGE " -lt $COV ]; then
		if echo "${EXPECTED_LONG_COVERAGE} < ${COV}" | bc -l | grep -q 1; then
			_polap_log1 "  No reduction of the long-read data because $EXPECTED_LONG_COVERAGE < $COV"
			_polap_log3 ln -s $(realpath "$LR") "$nfq_file"
			ln -s $(realpath "$LR") "$nfq_file"
		else
			_polap_log1 "SUGGESTION: you might want to increase the minimum read lengths because you have enough long-read data."
			local RATE=$(echo "scale=3; $COV/$EXPECTED_LONG_COVERAGE" | bc)
			# Compare value with 0
			if echo "${RATE} > 0" | bc -l | grep -q 1; then
				_polap_log1 "  sampling long-read data by $RATE ..."
				_polap_log1 "    $RATE <= target long-read genome coverage[$COV]/expected long-read genome coverage[$EXPECTED_LONG_COVERAGE] ..."
				local seed=${_arg_seed:-$RANDOM}
				_polap_log0 "  random seed for reducing the whole-genome assembly long-read data: ${seed}"
				# _polap_log3 "seqkit sample -p ${RATE} ${LR} -o ${nfq_file}"
				# seqkit sample -p "${RATE}" "${LR}" -o "${nfq_file}" >${_polap_output_dest} 2>&1
				_polap_log3_cmd "seqkit sample -p ${RATE} -s ${seed} ${LR} -o ${nfq_file} 2>${_polap_output_dest}"
				_polap_log2 "  ${nfq_file}: a reduced long-read data is created"
			else
				_polap_log0 "  target coverage: ${COV}"
				_polap_log0 "  long-read coverage: ${EXPECTED_LONG_COVERAGE}"
				_polap_log0 "  sampling rate is ${COV} / ${EXPECTED_LONG_COVERAGE} => ${RATE}"
				_polap_log0 "  genome size: ${EXPECTED_GENOME_SIZE}"
				_polap_log0 "  total long-read: ${LONG_TOTAL_LENGTH}"
				_polap_log0 "  Too large expected long-read coverage"
				_polap_log0 "  Expected genome size may be too small."
				die "ERROR: long-read sampling rate is not greater than 0."
			fi
		fi
	fi
	LR="${nfq_file}"

	# purge the long-read data of shorter than $MR bp
	check_file_existence "${LR}"
	_polap_log1 "  keeps long reads of length being at least $MR bp ..."
	_polap_log2 "  deletes $LRNK"
	_polap_log3 rm -f "$LRNK"
	rm -f "$LRNK"
	_polap_log3 seqkit seq --quiet -m "${MR}" --threads 4 "${LR}" -o "${LRNK}"
	seqkit seq --quiet -m "${MR}" --threads 4 "${LR}" -o "${LRNK}" >${_polap_output_dest} 2>&1
	_polap_log0 "  output: ${LRNK}"
	_polap_log2 "  deletes $nfq_file"
	_polap_log3 rm "$nfq_file"
	rm "$nfq_file"
	_polap_log0 "creating the statisics for the reduced long-read data ..."
	_polap_log3 seqkit stats -T "${LRNK}" ">${_polap_var_base_nk_fq_stats}"
	seqkit stats -T "${LRNK}" >"${_polap_var_base_nk_fq_stats}"
	_polap_log0 "  output: ${_polap_var_base_nk_fq_stats}"

	_polap_log1 "NEXT (for testing purpose only): $(basename "$0") flye1 -g 150000"
	_polap_log1 "NEXT (for testing purpose only): $(basename "$0") flye1 --test"
	_polap_log1 "NEXT: $(basename $0) flye1 -o $ODIR [-t $NT] [-c $COV]"

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}

################################################################################
# Executes Flye for a whole-genome assembly upto the contigger stage
#
# Arguments:
#   -t $NT: the number of CPU cores
#   -c $COV: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   $ODIR/short_expected_genome_size.txt (ignored with -g option)
#   $LRNK
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
	source "$script_dir/polap-variables-base.sh"
	source "$script_dir/polap-variables-wga.sh"
	source "$script_dir/run-polap-function-utilities.sh"

	LRNK="${_polap_var_base_nk_fq_gz}"

	help_message=$(
		cat <<HEREDOC
# Executes Flye for a whole-genome assembly upto the contigger stage
#
# Arguments:
#   -t $NT: the number of CPU cores
#   -c $COV: the Flye's coverage option
#   -g <arg>: computed by find-genome-size menu or given by users
# Inputs:
#   $ODIR/short_expected_genome_size.txt (ignored with -g option)
#   $LRNK
# Outputs:
#   $ODIR/0/30-contigger/contigs.fasta
#   $ODIR/0/30-contigger/contigs_stats.txt
#   $ODIR/0/30-contigger/graph_final.fasta
#   ${_polap_var_wga_contigger_gfa}
Example: $(basename $0) ${_arg_menu[0]} [-t|--threads <arg>] [-c|--coverage <arg>] [-g|--genomesize <arg>]
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" ]] && _polap_echo0 "${help_message}" && exit $EXIT_SUCCESS

	check_file_existence "$ODIR/short_expected_genome_size.txt"
	check_file_existence "${LRNK}"

	_polap_log0 "assembling the genome using the reduced long-read data ..."
	_polap_log0 "  input1: ${_polap_var_base_genome_size}"
	_polap_log0 "  input2: ${LRNK}"

	#   ${_polap_var_wga_contigger_gfa}
	if [ -s "${_polap_var_wga_contigger_gfa}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log0 "  found: ${_polap_var_wga_contigger_gfa}, so skipping the whole-genome assembly ..."
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
			_polap_log1 "OPTION: short reads expected genome size (bases) we use instead: ${_expected_genome_size}"
		fi

		_polap_log1 "  executing the whole-genome assembly using flye ... be patient!"
		_polap_log2 flye --nano-raw "${LRNK}" \
			--out-dir "${_polap_var_wga}" \
			--threads "${NT}" \
			--asm-coverage "${COV}" \
			--genome-size "${EXPECTED_GENOME_SIZE}" \
			--stop-after contigger
		flye --nano-raw "${LRNK}" \
			--out-dir "${_polap_var_wga}" \
			--threads "${NT}" \
			--asm-coverage "${COV}" \
			--genome-size "${EXPECTED_GENOME_SIZE}" \
			--stop-after contigger \
			>${_polap_output_dest} 2>&1
	fi

	_polap_log1 "  assembly graph in the flye contigger stage: ${_polap_var_wga_contigger_gfa}"
	_polap_log1 "NEXT: $(basename $0) blast-genome -o $ODIR [-i 0]"
	_polap_log1 "NEXT: $(basename $0) annotate -o $ODIR [-i 0]"

	_polap_log2 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$DEBUG" -eq 1 ] && set +x
}
