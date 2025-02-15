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
# Ensure that the current script is sourced only once
source "$script_dir/run-polap-function-include.sh"
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

# Function to extract the number from the filename
_polap_utility_extract_number_circular_path() {
	local filename="$1"
	# Extract the base name (remove leading paths)
	local basename=$(basename "$filename")
	# Extract the number using sed
	local number=$(echo "$basename" | sed -n 's/.*circular_path_\([0-9]*\)_concatenated.*/\1/p')
	echo "$number"
}

# Function to determine the minimum of two values in files.
# Each of the two files has a single line of a number.
# Call the function with the provided files
# short_file="short_total_length.txt"
# long_file="long_total_length.txt"
# minimum=$(_polap_minimum_two_values_in_files "$short_file" "$long_file")
#
# Output the result
# echo "The minimum value is: $minimum"
_polap_minimum_two_values_in_files() {
	# Read the first file and store the value
	local value1=$(<"$1")

	# Read the second file and store the value
	local value2=$(<"$2")

	# Compare the two values
	if ((value1 < value2)); then
		echo "$value1"
	else
		echo "$value2"
	fi
}

# Function to create a Bash array of _arg_disassemble_n values
# input1: the starting value
# input2: the ending value
# input3: the number of items in the result (default: 10)
# if the number of array is 1, then use the starting value only.
#
# Example usage
# _arg_disassemble_a=10000000
# _arg_disassemble_b=100000000
# _arg_disassemble_n=10
#
# result=($(_polap_array_numbers_between_two "$_arg_disassemble_a" "$_arg_disassemble_b" "$_arg_disassemble_n"))
# result=($(_polap_array_numbers_between_two "$_arg_disassemble_a" "$_arg_disassemble_b"))
# echo "Generated array: ${result[@]}"
_polap_array_numbers_between_two() {
	local start="$1"
	local end="$2"
	local count="${3:-10}"
	local array=()
	local step
	local i

	if ((start == end)); then
		# If start and end are the same, create an array of repeated numbers
		for ((i = 0; i < count; i++)); do
			array+=("$start")
		done
	elif ((count == 1)); then
		# If count is 1, create an array with a single item using the end value
		# array=("$start")
		array=("$end")
	else
		# Calculate the step size
		step=$(((end - start) / (count - 1)))

		# Populate the array with evenly spaced values
		for ((i = 0; i < count; i++)); do
			array+=($((start + i * step)))
		done

		# Ensure the last value is exactly _arg_disassemble_b
		array[$((count - 1))]="$end"
	fi

	echo "${array[@]}"
}

# Function to convert a comma-delimited string to a Bash array
# input="A_1.fastq,A_2.fastq,A_3.fastq"
# result=($(convert_elements_to_array "$input"))
_polap_convert_elements_to_array() {
	local input_string=$1
	IFS=',' read -r -a array <<<"$input_string"
	echo "${array[@]}"
}

# Function to handle FASTQ files
# output: fastq file
# input1: fastq file
# input2: fastq file (optional)
#
# Example usage
# _polap_concatenate_fastq_to "s1.fq" "file1.fq"
# _polap_concatenate_fastq_to "s1.fq" "file1.fq" "file2.fq.gz"
# _polap_concatenate_fastq_to "s1.fq" "file1.fq.gz" "file2.fq.gz"
_polap_concatenate_fastq_to() {
	local output_file="$1"     # Output file (e.g., s1.fq)
	local input_file1="$2"     # First input file
	local input_file2="${3:-}" # Second input file (optional)

	# Check if the first input file exists
	if [[ ! -f "$input_file1" ]]; then
		die "ERROR: no such file: $input_file1"
	fi

	# Create or overwrite the output file
	rm -f "${output_file}"

	# If only one file is provided and it is not compressed, create a soft link
	if [[ -z "$input_file2" && ! "$input_file1" == *.gz ]]; then
		ln -s "$(realpath "$input_file1")" "$output_file" || {
			die "ERROR: cannot create a soft link."
		}
		_polap_log2 "  soft link created: $output_file -> $input_file1"
		return 0
	fi

	# Function to process a single file
	(
		process_file() {
			local file="$1"
			if [[ $file == *.gz ]]; then
				gunzip -c "$file" >>"$output_file" || {
					die "ERROR: gunzip -c $file"
				}
			else
				cat "$file" >>"$output_file" || {
					die "ERROR: cat $file"
				}
			fi
		}

		# Process the first file
		process_file "$input_file1"

		# Process the second file if provided
		if [[ -n "$input_file2" ]]; then
			if [[ ! -f "$input_file2" ]]; then
				die "ERROR: no such file: $input_file2"
			fi
			process_file "$input_file2"
			_polap_log2 "  concatenated fastq file: $input_file1 + $input_file2 -> $output_file"
		else
			_polap_log2 "  concatenated fastq file: $input_file1 -> $output_file"
		fi
	)

}

# Function to adjust a value of _alpha
# Example usage
# _alpha=5.5
# _delta=2.0
#
# echo "Original _alpha: $_alpha"
#
# # Increase _alpha
# _alpha=$(_polap_adjust_alpha_by_delta "$_alpha" "$_delta" "increase")
# echo "Updated _alpha (after increase): $_alpha"
#
# # Decrease _alpha
# _alpha=$(_polap_adjust_alpha_by_delta "$_alpha" "$_delta" "decrease")
# echo "Updated _alpha (after decrease): $_alpha"
#
# # Decrease _alpha to become negative and test revert
# _alpha=1.0
# _delta=2.0
# _alpha=$(_polap_adjust_alpha_by_delta "$_alpha" "$_delta" "decrease")
# echo "Updated _alpha (after revert): $_alpha"
# Function to adjust a value of _alpha
_polap_adjust_alpha_by_delta() {
	local _alpha=$1     # Current value of _alpha (float)
	local _delta=$2     # Delta value (float)
	local _direction=$3 # Direction: "increase" or "decrease"

	if [[ "$_direction" == "increase" ]]; then
		# Increase _alpha by _delta
		_alpha=$(echo "scale=5; $_alpha + $_delta" | bc)
	elif [[ "$_direction" == "decrease" ]]; then
		# Decrease _alpha by _delta
		_alpha=$(echo "scale=5; $_alpha - $_delta" | bc)
		# Check if _alpha is negative
		if (($(echo "$_alpha < 0" | bc -l))); then
			_alpha=$(echo "scale=5; $_alpha + $_delta" | bc)
		fi
	else
		die "Invalid direction. Use 'increase' or 'decrease'."
		return 1
	fi

	echo "$_alpha" # Output the updated value
}

# count the number of bases of a short-read fastq file
# input1: a fastq file
# output: a text file with the lengeth
_polap_total-length-short() {
	local input1=$1
	local output=$2

	check_file_existence "${input1}"

	_polap_log3_pipe "seqkit stats -Ta ${input1} |\
		csvtk cut -t -f sum_len |\
		csvtk del-header \
		>${output}"

	return 0
}

# unzip a gzipped file leaving the input as is.
# do not unzip if the input is not a gzipped file.
#
# unzipped_file=$(_polap_gunzip_file "${_short_read1}")
# _rstatus="$?"
# if [[ "$_rstatus" -eq 0 ]]; then
#   _polap_log2 "  unzipped file: $unzipped_file"
#   _short_read1="$unzipped_file"
# fi
_polap_gunzip_file() {
	local input_file="$1"

	# Check if the input file exists
	if [[ ! -f "$input_file" ]]; then
		die "Error: File '$input_file' not found."
	fi

	# Check if the file is gzipped
	if file "$input_file" | grep -q "gzip compressed data"; then
		# Extract the file name without the .gz extension
		local output_file="${input_file%.gz}"

		# Unzip the file and keep the original
		if gunzip -c "$input_file" >"$output_file"; then
			echo "$output_file"
			return 0
		else
			die "ERROR: failed to unzip '$input_file'."
		fi
	else
		_polap_log2 "  file '$input_file' is not gzipped."
		return 1
	fi
}

# create index.txt (TSV)
# --disassemble-a
# --disassemble-b
# --disassemble-n
# --genomesize-a
# --genomesize-b
# --genomesize-n
# e.g., index sample-size genome-size
# 0	10000000	150000
# 1	14444444	150000
# 2	18888888	150000
# 3	23333332	150000
# 4	27777776	150000
# 5	32222220	150000
function _run_polap_step-disassemble-make-index-for-b {
	local _index_table
	_index_table="${1}"

	local index
	local i
	local j
	local _datasize_array
	local _datasize_array_length
	local sampling_datasize
	local genomesize_array
	local genomesize_array_length
	local genomesize

	# Add commands for step 4
	if [[ "${_arg_short_read1_is}" == "off" ]] &&
		[[ "${_arg_short_read2_is}" == "off" ]]; then
		local total_size_data=$(cat "${_polap_var_outdir_long_total_length}")
		_polap_log2 "  the minimum read data size is just the long-read (bp): $total_size_data"
	else
		local total_size_data=$(_polap_minimum_two_values_in_files \
			"${_polap_var_outdir_long_total_length}" \
			"${_polap_var_outdir_short_total_length}")
		_polap_log2 "  the minimum read data size of the long- and short-read (bp): $total_size_data"
	fi

	if [[ "${_arg_disassemble_b_is}" == "off" ]]; then
		_polap_log2 "  the largest read data size set to (bp): $total_size_data"
		_arg_disassemble_b=${total_size_data}
	else
		if ((total_size_data < _arg_disassemble_b)); then
			_polap_log2 "  the largest read data size set to (bp): $total_size_data"
			_arg_disassemble_b=${total_size_data}
		fi
	fi
	_polap_log1 "    the disassemble N: ${_arg_disassemble_n}"
	_polap_log1 "    the disassemble A (bp): ${_arg_disassemble_a}"
	_polap_log1 "    the disassemble B (bp): ${_arg_disassemble_b}"
	_polap_log1 "    the disassemble M (bp): ${_arg_disassemble_m}"
	_polap_log1 "    the disassemble alpha: ${_arg_disassemble_alpha}"
	_polap_log1 "    the disassemble delta: ${_arg_disassemble_delta}"
	if ((_arg_disassemble_a <= _arg_disassemble_b)); then
		_polap_log2 "  check: ${_arg_disassemble_a} is less than or equal to ${_arg_disassemble_b}."
	else
		die "ERROR: ${_arg_disassemble_a} (A) is greater than ${_arg_disassemble_b} (B)."
	fi

	_datasize_array=($(_polap_array_numbers_between_two \
		"$_arg_disassemble_a" \
		"$_arg_disassemble_b" \
		"$_arg_disassemble_n"))
	_polap_log2 "  generated step array (${#_datasize_array[@]}): ${_datasize_array[@]:0:3} ... ${_datasize_array[@]: -2}"

	# consider genomesize_array
	if [[ -n "${_arg_genomesize}" ]]; then
		_arg_genomesize_a="${_arg_genomesize}"
		_arg_genomesize_b="${_arg_genomesize}"
		_arg_genomesize_b_is="on"
		_arg_genomesize_n=1
	fi

	if [[ "${_arg_genomesize_b_is}" == "on" ]]; then
		genomesize_array=($(_polap_array_numbers_between_two \
			"$_arg_genomesize_a" \
			"$_arg_genomesize_b" \
			"$_arg_genomesize_n"))
		genomesize_array_length=${#genomesize_array[@]}
		_polap_log1 "  genome sizes (${_arg_genomesize_n}): ${_arg_genomesize_a} ~ ${_arg_genomesize_b}"
	fi

	_datasize_array_length=${#_datasize_array[@]}

	_polap_log3_cmd rm -f "${_index_table}"
	index=0
	for ((i = 0; i < _datasize_array_length; i++)); do
		sampling_datasize="${_datasize_array[$i]}"
		if [[ "${_arg_genomesize_b_is}" == "on" ]]; then
			for ((j = 0; j < genomesize_array_length; j++)); do
				genomesize="${genomesize_array[$j]}"
				printf "%d\t%d\t%d\n" $index $sampling_datasize $genomesize \
					>>"${_index_table}"
				index=$((index + 1))
			done
		else
			printf "%d\t%d\n" $i $sampling_datasize \
				>>"${_index_table}"
		fi
	done
	return 0
}

# create index.txt (TSV) using the proportion
# --disassemble-p
# --disassemble-n
# --genomesize-a
# --genomesize-b
# --genomesize-n
# e.g., index sample-size genome-size
# 0	10000000	150000
# 1	14444444	150000
# 2	18888888	150000
# 3	23333332	150000
# 4	27777776	150000
# 5	32222220	150000
function _run_polap_step-disassemble-make-index-for-p {
	local _index_table
	_index_table="${1}"

	local index
	local i
	local j
	local _datasize_array
	local _datasize_array_length
	local sampling_datasize
	local genomesize_array
	local genomesize_array_length
	local genomesize
	local disassemble_a
	local disassemble_b
	local total_size_data

	total_size_data=$(<"${_polap_var_outdir_long_total_length}")
	disassemble_b=$(_polap_utility_compute_percentage \
		"${_arg_disassemble_p}" \
		"$total_size_data")

	_polap_log1 "    the disassemble N: ${_arg_disassemble_n}"
	_polap_log1 "    the disassemble P (%): ${_arg_disassemble_p}"
	_polap_log1 "    the disassemble A (bp): ${_arg_disassemble_a}"
	_polap_log1 "    the disassemble B (bp): ${disassemble_b}"
	_polap_log1 "    the disassemble M (bp): ${_arg_disassemble_m}"
	_polap_log1 "    the disassemble alpha: ${_arg_disassemble_alpha}"
	_polap_log1 "    the disassemble delta: ${_arg_disassemble_delta}"
	if ((_arg_disassemble_a <= disassemble_b)); then
		_polap_log2 "  check: ${_arg_disassemble_a} is less than or equal to ${disassemble_b}."
	else
		die "ERROR: ${_arg_disassemble_a} (A) is greater than ${disassemble_b} (B)."
	fi

	_datasize_array=($(_polap_array_numbers_between_two \
		"$_arg_disassemble_a" \
		"$disassemble_b" \
		"$_arg_disassemble_n"))
	_polap_log2 "  generated step array (${#_datasize_array[@]}): ${_datasize_array[@]:0:3} ... ${_datasize_array[@]: -2}"

	# consider genomesize_array
	if [[ -n "${_arg_genomesize}" ]]; then
		_arg_genomesize_a="${_arg_genomesize}"
		_arg_genomesize_b="${_arg_genomesize}"
		_arg_genomesize_b_is="on"
		_arg_genomesize_n=1
	fi

	if [[ "${_arg_genomesize_b_is}" == "on" ]]; then
		genomesize_array=($(_polap_array_numbers_between_two \
			"$_arg_genomesize_a" \
			"$_arg_genomesize_b" \
			"$_arg_genomesize_n"))
		genomesize_array_length=${#genomesize_array[@]}
		_polap_log1 "  genome sizes (${_arg_genomesize_n}): ${_arg_genomesize_a} ~ ${_arg_genomesize_b}"
	fi

	_datasize_array_length=${#_datasize_array[@]}

	_polap_log3_cmd rm -f "${_index_table}"
	index=0
	for ((i = 0; i < _datasize_array_length; i++)); do
		sampling_datasize="${_datasize_array[$i]}"
		if [[ "${_arg_genomesize_b_is}" == "on" ]]; then
			for ((j = 0; j < genomesize_array_length; j++)); do
				genomesize="${genomesize_array[$j]}"
				printf "%d\t%d\t%d\n" $index $sampling_datasize $genomesize \
					>>"${_index_table}"
				index=$((index + 1))
			done
		else
			printf "%d\t%d\n" $i $sampling_datasize \
				>>"${_index_table}"
		fi
	done
	return 0
}

# Estimate the genome size from a short-read dataset
#
# input1: a short-read fastq data
# input2: an output folder _outdir
# input3: redo default is on
# output: the genome size
# _outdir_genome_size="${_outdir}/short_expected_genome_size.txt"
_polap_find-genome-size() {
	local _short_read1="$1"
	local _outdir="$2"
	local _redo="${3:-on}"
	local _rstatus
	local _outdir_genome_size
	local _outdir_jellyfish_out
	local _outdir_jellyfish_out_histo
	local unzipped_file
	local _EXPECTED_GENOME_SIZE
	local _expected_genome_size_bp

	# Set paths for bioproject data
	# source "$script_dir/polap-variables-common.sh" # '.' means 'source'
	# source "$script_dir/run-polap-function-utilities.sh"

	_polap_log1 "  estimate the genome size using the short-read data ${_short_read1} ..."

	if [[ -d "${_outdir}" ]]; then
		_polap_log2 "  output folder: ${_outdir}"
	else
		_polap_log3 mkdir -p "${_outdir}"
		mkdir -p "${_outdir}"
	fi
	check_file_existence "${_short_read1}"

	_polap_log2 "  input1: ${_short_read1}"
	_outdir_genome_size="${_outdir}/short_expected_genome_size.txt"
	_outdir_jellyfish_out="${_outdir}/jellyfish_out"
	_outdir_jellyfish_out_histo="${_outdir}/jellyfish_out.histo"

	# See https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/
	if [ -s "${_outdir_genome_size}" ] && [ "${_redo}" = "off" ]; then
		_polap_log2 "  found: ${_outdir_genome_size}, so skipping the genome size estimation ..."
		_polap_log3_file "${_outdir_genome_size}"
	else
		if [ -s "${_outdir_jellyfish_out}" ] && [ "${_redo}" = "off" ]; then
			_polap_log2 "  found: ${_outdir_jellyfish_out}, so skipping the JellyFish counting ..."
			_polap_log3_file "${_outdir_jellyfish_out}"
		else
			unzipped_file=$(_polap_gunzip_file "${_short_read1}")
			_rstatus="$?"
			if [[ "$_rstatus" -eq 0 ]]; then
				_polap_log2 "  unzipped file: $unzipped_file"
				_short_read1="$unzipped_file"
			fi

			if [[ -s "${_short_read1}" ]]; then
				_polap_log3_pipe "jellyfish count \
					-t ${_arg_threads} -C -m 19 \
					-s 5G \
					-o ${_outdir_jellyfish_out} \
					--min-qual-char=? \
					${_short_read1}"
			else
				die "ASSERT: we must have at least one short-read fastq file."
			fi
			check_file_existence "${_outdir_jellyfish_out}"
		fi

		check_file_existence "${_outdir_jellyfish_out}"
		_polap_log3_pipe "jellyfish histo \
			-o ${_outdir_jellyfish_out_histo} \
			${_outdir_jellyfish_out}"
		_polap_log3_file "${_outdir_jellyfish_out_histo}"

		_polap_log3_pipe "Rscript --vanilla $script_dir/run-polap-r-jellyfish.R \
			${_outdir_jellyfish_out_histo} \
			${_outdir_genome_size}"
		# Check the exit status
		if [ $? -ne 0 ]; then
			# Take action if needed, e.g., logging, sending a notification, etc.
			die "ERROR: not being to estimate the genome size using short-read data: short-read may be too small."
		fi
	fi
	_polap_log2 "  output: ${_outdir_genome_size}"
	_polap_log3_cat "${_outdir_genome_size}"

	local _EXPECTED_GENOME_SIZE=$(<"${_outdir_genome_size}")
	local _EXPECTED_GENOME_SIZE=${_EXPECTED_GENOME_SIZE%.*}
	local _expected_genome_size_bp=$(_polap_utility_convert_bp ${_EXPECTED_GENOME_SIZE})
	_polap_log2 "  expected genome size using short-read data (bases): ${_expected_genome_size_bp}"

	return 0
}

# src/polap.sh disassemble-stats -o jvalidus2 --anotherdir jvalidus
function _run_polap_disassemble-stats {
	source "$script_dir/polap-variables-common.sh"
	source "$script_dir/polap-function-disassemble-seeds.sh"

	help_message=$(
		cat <<HEREDOC
Show some stats of disassemble menu results.
HEREDOC
	)

	_arg_plastid="on"
	local _input_short_reads="${_arg_outdir}/s1.fq"
	local _disassemble_dir="${_arg_outdir}/disassemble"

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return 0
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	local distribution1_file="${_disassemble_dir}/summary.txt"
	local distribution2_file="${_disassemble_dir}/x/summary1.txt"
	local percentile_file="${_disassemble_dir}/summary-length-percentile.txt"
	_polap_log3_pipe "Rscript $script_dir/run-polap-r-disassemble-stats.R \
			--dista ${distribution1_file} \
			--distb ${distribution2_file} \
			--out ${percentile_file}"

	# Check the exit status
	if [ $? -ne 0 ]; then
		# Take action if needed, e.g., logging, sending a notification, etc.
		die "ERROR: compare lengths in ${distribution1_file} vs. ${distribution2_file}"
	else
		_polap_log0_cat "${percentile_file}"
	fi

	return 0
}

# _run_polap_step-disassemble-cflye <out> <long.fq> <genomesize> <alpha> <resume:on>
#
function _run_polap_step-disassemble-cflye {
	local out="${1}"
	local long_read="${2}"
	local expected_genome_size="${3}"
	local alpha="${4}"
	local resume="${5:-off}"
	local _command1

	if [[ "${_POLAP_RELEASE}" -eq 0 ]]; then
		_command1="/home/goshng/all/polap/Flye/bin/cflye"
	else
		_command1="cflye"
	fi

	_command1+=" \
          ${_arg_flye_data_type} \
          ${long_read} \
			    --out-dir ${out} \
			    --disjointig-min-coverage ${alpha} \
			    --threads ${_arg_threads}"

	_command1+=" \
			  --asm-coverage ${_arg_flye_asm_coverage} \
			  --genome-size ${expected_genome_size}"
	if [[ "${_arg_flye_asm_coverage}" -le 0 ]]; then
		die "ERROR: must be positive --flye-asm-coverage ${_arg_flye_asm_coverage}"
	fi

	if [[ "${_arg_contigger}" == "on" ]]; then
		_command1+=" \
		    --stop-after contigger"
	fi

	if [[ "${resume}" != "off" ]]; then
		_command1+=" \
		    --resume"
	fi

	if [[ "${_arg_debug}" == "on" ]]; then
		_command1+=" \
		      --debug"
	fi

	_command1+=" \
		    2>${_polap_output_dest}"

	_polap_log3_pipe "${_command1}"

	return 0
}

function _run_polap_step-disassemble-archive-cfile {
	local cfile="${1}"

	if [[ -s "${cfile}" ]]; then
		_polap_log3_cmd mkdir -p ${_arg_archive}/$(dirname "${cfile#*/}")
		_polap_log3_pipe "cp -p \
      ${cfile} \
      ${_arg_archive}/${cfile#*/}"
	fi
	return 0
}

################################################################################
# It is not as straightforward as it might seem. The simplest approach is to test plastid genome assembly by subsampling long-read data. A more complex approach involves using a slightly modified version of Flye, called cflye, which provides coverage options for generating disjointigs. Short-read data is still used to get a rough estimate of the overall genome size, but once subsampling is applied, the size of the genome assembled from long reads can change. The sizes of both the long-read and short-read data are defined by their total number of bases. When subsampling, we reduce the amount of short-read data proportionally to match the size of the subsampled long-read data, then use JellyFish (a k-mer counter) to predict the total size of the genome assembled from these long reads.
#
# If the total long-read data is smaller than the short-read data, we simply subsample the short reads to match the long-read dataset size. However, if the long-read data is larger than the short-read data, we use all the short reads. Because our goal is to estimate the approximate genome size that can be assembled from subsampled long reads, some level of error is to be expected.
#
# In the Flye genome assembly program, disjointigs are generated, and those whose read coverage falls below a certain cutoff are filtered out. We set this initial cutoff using the size of the long-read data (A = 30 Mb, where Mb stands for megabases) and alpha (α = 0.1). We then compare the total length (S) of the disjointigs against a fixed value M. If S is greater than M, alpha is increased by an amount delta (δ = 0.1); if S is less than M, alpha is decreased by the same delta. Along with this process, after each cycle we also increase the long-read data size A by d (for example, 10 Mb) and repeat these cycles until A reaches B, for a total of n cycles. The values of A, B, and n are predetermined, and d is derived from these values. For instance, if A = 30 Mb, B = 120 Mb, and n = 10, then d = 10 Mb.
#
# At the end of each cycle, we extract and evaluate potential plastid genome candidates from the newly assembled genome. To qualify as a plastid candidate, the assembly must consist of a single circular contig or three contigs that form a circular structure. We then assess the size of this circular contig, the number of expected genes, and so on. These values can be plotted as Y-coordinates against the cycle number on the X-axis in a scatter plot, allowing us to visualize progress across cycles.
#
#### Step 1: Estimating Long-Read Dataset Size
#
# The process starts by determining the total length of the long-read dataset by counting the total number of bases present in the sequencing reads. This information is crucial as it defines the upper limit for subsampling and subsequent analyses.
#
#### Step 2 & 3: Preparing Short-Read Data or presteps 1 and 2
#
# When short-read data is available, the sequencing reads from different files are first combined into a single dataset to ensure uniformity. The total number of bases in this combined short-read dataset is then calculated. This base count is used to balance the amount of short-read data with the subsampled long-read data in later steps.
#
#### Step 4: Defining Subsampling Range
#
# A range of data sizes is defined to ensure a comprehensive assembly across varying input sizes. The minimum dataset size is set by comparing the total base counts of the long-read and short-read datasets, selecting the smaller of the two. If a maximum subsample size is not provided, it defaults to the size of the smaller dataset. This range is then divided into evenly spaced intervals, creating a systematic plan for subsampling and subsequent assembly.
#
#### Step 5: Subsampling Long-Read Data
#
# To explore the effect of varying coverage levels on assembly quality, smaller subsets of the long-read dataset are created. The proportion of reads to be included in each subset is determined by dividing the target size by the total dataset size. Random selection of reads ensures unbiased sampling, and the process is repeated across the predefined range of data sizes.
#
#### Step 6: Subsampling Short-Read Data
#
# Short-read data is similarly subsampled to match the size of the corresponding long-read subset. This ensures balanced input for assembly, which is important for accurate genome size estimation and assembly quality.
#
#### Step 7: Genome Size Estimation
#
# If a reference genome size is not available, the genome size is estimated using the subsampled short-read data. A k-mer counting tool is employed to generate a histogram of k-mer frequencies, which is then analyzed to derive an estimated genome size. This value serves as a reference for determining coverage thresholds in subsequent assembly steps.
#
#### Step 8: Assembly Using Modified Genome Assembler
#
# The assembly process uses a specialized version of a genome assembler tailored for handling plastid genomes. Parameters such as the estimated genome size and minimum coverage thresholds are supplied to guide the assembly. The assembler produces intermediate outputs known as disjointigs, representing contiguous sequences before final polishing. Memory usage is monitored to prevent excessive resource consumption.
#
#### Step 9: Adjusting Coverage Threshold
#
# Post-assembly, the total length of the assembled genome is compared to a predefined threshold. If the length exceeds this threshold, the coverage requirement is increased to filter out low-coverage regions in the next iteration. Conversely, if the length is below the threshold, the coverage requirement is reduced. This iterative adjustment helps refine the assembly.
#
#### Step 10: Generating Assembly Graph Statistics
#
# Key statistics of the assembly graph, such as the number of contigs and total segment length, are calculated. These metrics provide an overview of the structural complexity of the assembly and help in evaluating its completeness.
#
#### Step 11: Annotating Contigs
#
# The assembled contigs are examined for the presence of plastid-specific genes. Contigs with a high density of plastid genes are prioritized for further analysis, while those containing predominantly non-plastid genes are excluded.
#
#### Step 12: Selecting Seed Contigs
#
# Seed contigs are selected based on specific criteria to identify candidates most likely to represent the plastid genome. This process involves several substeps:
#
# - **Step 12-1: Preselection of Contigs**
#   Contigs are initially filtered based on gene density. Contigs with a sufficient number of plastid genes and minimal contamination from non-plastid genes are retained.
#
# - **Step 12-2: Determining Depth Range**
#   The sequencing depth of selected contigs is analyzed, and an appropriate depth range is defined. Contigs outside this range are excluded to ensure that only those with typical plastid genome coverage are retained.
#
# - **Step 12-3: Filtering the Assembly Graph**
#   The assembly graph is filtered by removing contigs and connections that fall outside the selected depth range. This step simplifies the graph and focuses on high-confidence regions.
#
# - **Step 12-4: Identifying Connected Components**
#   The filtered graph is analyzed to find clusters of interconnected contigs. These clusters represent potential plastid genome sequences. The contigs forming these clusters are selected as seed contigs for further analysis.
#
#### Step 13: Extracting Circular Genome Sequences
#
# Circular genome sequences are identified by tracing paths through the assembly graph that form closed loops. Only complete circular structures are retained, as plastid genomes are typically circular. The extracted sequences are saved for further evaluation.
#
#### Step 14: Selecting the Best Genome Candidate
#
# When a reference plastid genome is available, the extracted candidates are compared against it. The candidate with the highest sequence similarity and coverage is selected as the final assembly. In cases where no reference is provided, the first complete circular sequence is chosen.
#
#### Step 15: Finalizing the Assembly
#
# The final plastid genome assembly is validated by calculating key metrics such as total length and GC content. These metrics are logged, and the assembled sequence is saved as the final result. This completes the iterative workflow, yielding a high-quality plastid genome assembly.
#
################################################################################
function _run_polap_disassemble {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_/Menu_/)"

	local i
	local j
	local k
	local c
	local n
	local s
	local d
	local col1
	local cfile
	local _rstatus
	local index
	local rate
	local size
	local alpha
	local genomesize
	local output
	local long_read
	local randomseed
	local _disassemble_dir
	local _summary_table
	local _input_short_reads
	local _ptdna
	local _ptdir
	local _ref_ptdna
	local summary_j_ordered
	local _disassemble_i
	local summary_table
	local _index_table
	local _mtcontigname
	local _include
	local _exclude
	local _stage_array
	local _disassemble_i
	local _disassemble_i_stage
	local _summary1_ordered
	local _disassemble_params_file

	source "$script_dir/polap-variables-common.sh"
	local _disassemble_dir="${_arg_outdir}/disassemble"
	local _summary_table="${_disassemble_dir}/summary.txt"

	_arg_plastid="on"
	_arg_disassemble_s=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_s}")
	_arg_disassemble_a=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_a}")
	_arg_disassemble_b=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_b}")
	_arg_disassemble_m=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_m}")
	_arg_disassemble_s_max=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_s_max}")
	_arg_genomesize_a=$(_polap_utility_convert_unit_to_bp "${_arg_genomesize_a}")
	_arg_genomesize_b=$(_polap_utility_convert_unit_to_bp "${_arg_genomesize_b}")

	if [[ -n "${_arg_genomesize}" ]]; then
		_arg_genomesize=$(_polap_utility_convert_unit_to_bp "${_arg_genomesize}")
	fi
	local _input_short_reads="${_arg_outdir}/s1.fq"

	help_message=$(
		cat <<HEREDOC
Plastid genome assembly by subsampling long-read data without references

Inputs
------

- long-read data: ${_arg_long_reads} (default: l.fq)
- short-read data 1 : ${_arg_short_read1} (no default)
- short-read data 2 (optional): ${_arg_short_read2} (no default)

--disassemble-p ${_arg_disassemble_p}: the percentile of the largest long read
--disassemble-s: the threshold to change alpha
--disassemble-stop-after: stage1, stage2, stage3

Outputs
-------

- plastid genome assembly: ${_arg_final_assembly}
- trace plots for the features of plastid genome assemblies
  ${_arg_outdir}/disassemble/0/summary1-ordered.pdf
  ${_arg_outdir}/disassemble/0/summary1-ordered.txt
  ${_arg_outdir}/disassemble/x/summary1-ordered.pdf
  ${_arg_outdir}/disassemble/x/summary1-ordered.txt

Arguments
---------

Stage 1:
-o ${_arg_outdir}
-l ${_arg_long_reads}: a long-read fastq data file
-a ${_arg_short_read1}: a short-read fastq data file 1
-b ${_arg_short_read2} (optional): a short-read fastq data file 2
--disassemble-min-memory ${_arg_disassemble_min_memory}: the minimum memory in Gb
--disassemble-a ${_arg_disassemble_a}: the smallest long read
--disassemble-b ${_arg_disassemble_b}: the largest long read
--disassemble-n ${_arg_disassemble_n}: the number of steps (max 1000)
--disassemble-i ${_arg_disassemble_i}: the index of the stage 1
--disassemble-p ${_arg_disassemble_p}: the percentile of the largest long read
--disassemble-m ${_arg_disassemble_m}: the threshold to change alpha
--disassemble-alpha ${_arg_disassemble_alpha}: the minimum Flye's disjointig coverage
--disassemble-delta ${_arg_disassemble_delta}: the move size of alpha (0.1 - 1.0)
--disassemble-compare-to-fasta <FASTA>
--start-index <INDEX>: from <INDEX>
--end-index <INDEX>: to <INDEX> - 1
--steps-include <STPES>: STEPS can be 1,2,3 or 1-15
--steps-exclude <STPES>: STEPS can be 1,2,3 or 1-15
-t ${_arg_threads}: the number of CPU cores
--random-seed <arg>: 5-digit number
--no-contigger for a single final step: use the long-read polished gfa


--stages-include <STAGES> for internal use
--disassemble-stop-after stage1 for external interface
--disassemble-stop-after stage2 for external interface

--genomesize <NUMBER>: the number of steps is equal to 1.
--genomesize-a ${_arg_genomesize_a}: the smallest genome-size step size
--genomesize-b ${_arg_genomesize_b}: the largest genome-size step size
--genomesize-n ${_arg_genomesize_n}: the number of steps in the genome-size



Stage 2:
--disassemble-s: the threshold to change alpha
--disassemble-alpha ${_arg_disassemble_alpha}: the minimum Flye's disjointig coverage
--disassemble-r ${_arg_disassemble_r}: the number of replicates

Test stage: with neither of these two options, one sequence is selected as 
a final genome.
--species <SPECIES>
--disassemble-compare-to-fasta <FASTA>

Menus
-----

- help: display this help message
- redo: overwrite previous results
- view: show some results
- reset or clean: delete output files 

Usages
------
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads}
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads} -a ${_arg_short_read1}
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads} -a ${_arg_short_read1} -b ${_arg_short_read2}

Examples
--------
1:
$(basename "$0") x-ncbi-fetch-sra --sra SRR7153095
$(basename "$0") x-ncbi-fetch-sra --sra SRR7161123
$(basename "$0") get-mtdna --plastid --species "Eucalyptus pauciflora"
cp o/00-bioproject/2-mtdna.fasta ptdna-Eucalyptus_pauciflora-known.fa

bash ptgaul/ptGAUL.sh -o o-ptgaul -r ptdna-Eucalyptus_pauciflora-known.fa -g 180000 -l long.fastq
cp -pr o-ptgaul/result_3000 o

2:
Stage 1 - long-read and short-read data
$(basename "$0") disassemble -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq --disassemble-i 1 --disassemble-p 1 --disassemble-n 10

$(basename "$0") disassemble archive

$(basename "$0") disassemble -o o-a --disassemble-best --disassemble-compare-to-fasta ptdna-epauciflora.fa --disassemble-i 1 -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq 

$(basename "$0") disassemble view
3:
$(basename "$0") disassemble view 2
$(basename "$0") disassemble -o o2 --anotherdir o --disassemble-s 314588465 --disassemble-alpha 3.5 -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq

For ptDNA of Juncus validus
$(basename "$0") get-mtdna -o jvalidus --plastid --species "Juncus validus"
$(basename "$0") disassemble -o jvalidus -l SRR21976089.fastq -a SRR21976091_1.fastq -b SRR21976091_2.fastq --disassemble-min-memory 9 -v --disassemble-a 50mb --disassemble-n 20
src/polap.sh disassemble -l SRR21976089.fastq -a SRR21976091_1.fastq -b SRR21976091_2.fastq -o jvalidus --disassemble-min-memory 9 -v --disassemble-a 50mb --disassemble-n 20

Stage 1 - long-read only
$(basename "$0") disassemble -l SRR7153095.fastq --disassemble-compare-to-fasta ptdna-epauciflora.fa --disassemble-min-memory 9

Menu examples:
$(basename "$0") disassemble view
$(basename "$0") disassemble view 0
$(basename "$0") disassemble view x
$(basename "$0") disassemble view best 0
$(basename "$0") disassemble view best x
$(basename "$0") disassemble report
$(basename "$0") disassemble report coverage <- with known ptDNA
$(basename "$0") disassemble best 0 46
$(basename "$0") disassemble best x 2
$(basename "$0") disassemble archive
$(basename "$0") disassemble archive polishing <- backup msbwt as well
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return 0
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		case "${_arg_menu[2]}" in
		best)
			_polap_log0_column "${_disassemble_dir}/${_arg_menu[3]}/summary1-ordered.txt"
			;;
		x)
			_polap_log0_column "${_disassemble_dir}/${_arg_menu[2]}/summary1.txt"
			;;
		'' | *[!0-9]*)
			_polap_log0_column "${_disassemble_dir}/0/summary1.txt"
			;;
		*)
			_polap_log0_column "${_disassemble_dir}/${_arg_menu[2]}/summary1.txt"
			;;
		esac

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# menu: report
	if [[ "${_arg_menu[1]}" == "report" ]]; then

		# report 1 --disassemble-i 3
		if [[ "${_arg_menu[2]}" == "1" ]]; then
			s=${_disassemble_dir}/${_arg_disassemble_i}/1/summary1.txt
			d=${_disassemble_dir}/${_arg_disassemble_i}/1/summary1.md
			csvtk -t cut -f index,long_rate_sample,alpha,peak_ram_size_gb,expected_genome_size,time,gfa_number_segments,gfa_total_segment_length,num_circular_paths,coverage_ref,coverage_target,length \
				"${s}" |
				csvtk -t round -n 4 -f long_rate_sample,coverage_ref,coverage_target |
				csvtk -t round -n 2 -f alpha |
				csvtk rename -t -f index,long_rate_sample,alpha,peak_ram_size_gb,expected_genome_size,time,gfa_number_segments,gfa_total_segment_length,num_circular_paths,coverage_ref,coverage_target,length \
					-n I,Rate,Alpha,'Pmem','G',Time,'N','L','C','C ref','C target',Length |
				csvtk -t csv2md -a right - \
					>"${d}"
			_polap_log0_head "${d}"
			d1="${_disassemble_dir}/${_arg_disassemble_i}/1/summary1-ordered.txt"
			d2="${_disassemble_dir}/${_arg_disassemble_i}/1/summary1-ordered.pdf"
			_polap_log3_pipe "Rscript --vanilla $script_dir/run-polap-r-disassemble.R \
        --table ${s} \
        --out ${d1} \
        --plot ${d2}"
		fi

		# report 2 --disassemble-i 3
		if [[ "${_arg_menu[2]}" == "2" ]]; then
			# concatenate all summary1.txt to the summary table.
			s=${_disassemble_dir}/${_arg_disassemble_i}/2/summary1.txt
			d=${_disassemble_dir}/${_arg_disassemble_i}/2/summary1.md
			csvtk -t cut -f index,long_rate_sample,alpha,peak_ram_size_gb,expected_genome_size,time,gfa_number_segments,gfa_total_segment_length,num_circular_paths,coverage_ref,coverage_target,length,pident \
				"${s}" |
				csvtk -t round -n 4 -f long_rate_sample,coverage_ref,coverage_target |
				csvtk -t round -n 2 -f alpha |
				csvtk rename -t -f index,long_rate_sample,alpha,peak_ram_size_gb,expected_genome_size,time,gfa_number_segments,gfa_total_segment_length,num_circular_paths,coverage_ref,coverage_target,length,pident \
					-n I,Rate,Alpha,'Pmem','G',Time,'N','L','C','C ref','C target',Length,'Identity' |
				csvtk -t csv2md -a right - \
					>"${d}"
			_polap_log0_head "${d}"
			d1="${_disassemble_dir}/${_arg_disassemble_i}/2/summary1-ordered.txt"
			d2="${_disassemble_dir}/${_arg_disassemble_i}/2/summary1-ordered.pdf"
			_polap_log3_pipe "Rscript --vanilla $script_dir/run-polap-r-disassemble.R \
        --table ${s} \
        --out ${d1} \
        --plot ${d2}"
		fi

		# report 3 --disassemble-i 3
		if [[ "${_arg_menu[2]}" == "3" ]]; then
			# concatenate all summary1.txt to the summary table.
			s0="${_disassemble_dir}/${_arg_disassemble_i}/1/summary1.txt"
			d1="${_disassemble_dir}/${_arg_disassemble_i}/1/summary1-scatter.txt"
			d2="${_disassemble_dir}/${_arg_disassemble_i}/1/summary1-scatter.pdf"
			_polap_log3_pipe "Rscript --vanilla $script_dir/run-polap-r-disassemble.R \
        --table ${s0} \
        --out ${d1} \
        --plot ${d2} \
			  --scatter"
		fi

		if [[ "${_arg_menu[2]}" == "9" ]]; then
			# concatenate all summary1.txt to the summary table.
			for i in "${_disassemble_dir}"/*/; do
				if [ -d "$i" ]; then
					i=$(basename "${i%/}")
					if [[ -s "${_disassemble_dir}/${i}/summary1.txt" ]]; then
						if [[ "${_arg_menu[2]}" == "coverage" ]]; then
							_polap_log3_pipe "Rscript --vanilla $script_dir/run-polap-r-disassemble.R \
        --table ${_disassemble_dir}/${i}/summary1.txt \
        --out ${_disassemble_dir}/${i}/summary1-ordered.txt \
        --plot ${_disassemble_dir}/${i}/summary1-ordered.pdf \
			  --coverage"
						else
							_polap_log3_pipe "Rscript --vanilla $script_dir/run-polap-r-disassemble.R \
        --table ${_disassemble_dir}/${i}/summary1.txt \
        --out ${_disassemble_dir}/${i}/summary1-ordered.txt \
        --plot ${_disassemble_dir}/${i}/summary1-ordered.pdf"
						fi
						_polap_log0_column ${_disassemble_dir}/${i}/summary1-ordered.txt
					fi
				fi
			done
		fi

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# menu: best
	if [[ "${_arg_menu[1]}" == "best" ]]; then
		if [[ "${_arg_menu[2]}" == "outfile" ]]; then
			_disassemble_i="${_disassemble_dir}/x"
		else
			_disassemble_i="${_disassemble_dir}/${_arg_menu[2]}"
		fi

		local j="${_arg_menu[3]}"
		if [[ "${_arg_menu[3]}" == "thirdfile" ]]; then
			local _summary1_ordered="${_disassemble_i}/summary1-ordered.txt"
			if [[ -s "${_summary1_ordered}" ]]; then
				j=$(awk 'NR==2 {print $1}' "${_summary1_ordered}")
			else
				_polap_log0 "ERROR: no such file: ${_summary1_ordered}"
			fi
		fi

		_ref_ptdna="${_disassemble_i}/${j}/ptdna.fa"

		summary_j_ordered="${_disassemble_i}/summary-${j}-ordered.txt"
		printf "%s\t%s\n" \
			"index" \
			"coverage" \
			>"${summary_j_ordered}"
		if [[ -s "${_ref_ptdna}" ]]; then
			for ((i = 0; i < 100; i++)); do
				local _ptdna="${_disassemble_i}/${i}/ptdna.fa"
				local _ptdir="${_disassemble_i}/${j}/best/${i}"
				mkdir -p "${_ptdir}"
				if [[ -s "${_ptdna}" ]]; then
					_polap_log3_pipe "python $script_dir/run-polap-py-compare2ptdna.py \
						    --seq1 ${_ref_ptdna} \
						    --seq2 ${_ptdna} \
					     --out ${_ptdir} \
					     2>$_polap_output_dest"
					c=$(<"${_ptdir}/coverage.txt")
				else
					_polap_log2 "  no such ptDNA assembled: ${_ptdna}"
					c=0
				fi
				printf "%s\t%s\n" \
					"${i}" \
					"${c}" \
					>>"${summary_j_ordered}"
			done
		else
			_polap_log0 "ERROR: no such file: ${_ref_ptdna}"
		fi

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# menu: archive
	# archive the disassemble analysis
	if [[ "${_arg_menu[1]}" == "archive" ]]; then
		_polap_log0 "archive ${_arg_outdir} to ${_arg_archive}"

		# base directory
		_polap_log3_cmd mkdir -p ${_arg_archive}
		_polap_log3_cmd cp -p ${_arg_outdir}/${_arg_log} ${_arg_archive}/${_arg_log}

		s="${_polap_var_outdir_long_total_length}"
		d="${_arg_archive}/${_polap_var_outdir_long_total_length#*/}"
		_polap_log3_cmd cp -p "${s}" "${d}"

		s="${_polap_var_outdir_short_total_length}"
		d="${_arg_archive}/${_polap_var_outdir_short_total_length#*/}"
		_polap_log3_cmd cp -p "${s}" "${d}"

		# 00-bioproject
		_polap_log3_cmd mkdir -p ${_arg_archive}/${_polap_var_project#*/}
		s="${_polap_var_project_mtdna_fasta1_stats}"
		d="${_arg_archive}/${_polap_var_project_mtdna_fasta1_stats#*/}"
		_polap_log3_cmd cp -p "${s}" "${d}"

		s="${_polap_var_project_mtdna_fasta2_accession}"
		d="${_arg_archive}/${_polap_var_project_mtdna_fasta2_accession#*/}"
		_polap_log3_cmd cp -p "${s}" "${d}"

		# msbwt
		if [[ "${_arg_menu[2]}" == "polishing" ]]; then
			_polap_log3_pipe "cp -pr \
        ${_polap_var_outdir_msbwt_dir} \
        ${_arg_archive}/${_polap_var_outdir_msbwt_dir#*/}"
		fi

		# disassemble
		_polap_log3_cmd mkdir -p ${_arg_archive}/${_disassemble_dir#*/}
		cp -p "${_arg_outdir}/timing-ptgaul.txt" "${_arg_archive}/"
		cp -p "${_arg_outdir}/ptdna-ptgaul.fa" "${_arg_archive}/"
		cp -p "${_arg_outdir}/ptdna-reference.fa" "${_arg_archive}/"

		# disassemble/0
		for i in "${_disassemble_dir}"/*/; do
			if [ -d "$i" ]; then
				i=$(basename "${i%/}")
				cp -p "${_arg_outdir}/timing-${i}.txt" "${_arg_archive}/"
				_disassemble_i="${_arg_outdir}/disassemble/${i}"

				_polap_log3_cmd mkdir -p ${_disassemble_i}
				mkdir -p "${_arg_archive}/disassemble/${i}"

				s="${_disassemble_i}/pt.1.fa"
				d="${_arg_archive}/${s#*/}"
				_polap_log3_cmd cp -p "${s}" "${d}"

				s="${_disassemble_i}/params.txt"
				d="${_arg_archive}/${s#*/}"
				_polap_log3_cmd cp -p "${s}" "${d}"

				for j in 1 2; do
					k="${_disassemble_dir}/${i}/${j}"
					mkdir -p "${_arg_archive}/disassemble/${i}/${j}"

					summary_table="${k}/summary1.txt"
					_index_table="${k}/index.txt"

					s="${k}"
					d="${_arg_archive}/${k#*/}"
					_polap_log3_cmd cp -p "${s}/summary"* "${d}/"
					_polap_log3_cmd cp -p "${s}/index.txt" "${d}/"

					# index.txt
					# Read the first column and iterate over it
					while IFS=$'\t' read -r col1 _; do
						# Add your commands to process $col1 here
						cfile="${k}/${col1}/30-contigger/graph_final.gfa"
						_run_polap_step-disassemble-archive-cfile "${cfile}"
						cfile="${k}/${col1}/30-contigger/graph_final.fasta"
						_run_polap_step-disassemble-archive-cfile "${cfile}"
						cfile="${k}/${col1}/mt.contig.name"
						_run_polap_step-disassemble-archive-cfile "${cfile}"
						cfile="${k}/${col1}/assembly_graph.gfa"
						_run_polap_step-disassemble-archive-cfile "${cfile}"
						cfile="${k}/${col1}/cflye.log"
						_run_polap_step-disassemble-archive-cfile "${cfile}"
						cfile="${k}/${col1}/params.json"
						_run_polap_step-disassemble-archive-cfile "${cfile}"
					done <"$_index_table"
				done

			fi
		done

		# ptdna
		# _polap_log3_pipe "cp ${_arg_outdir}/pt.1.fa \
		#     ${_arg_archive}/pt.1.fa"

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# menu: ptgaul
	# ptgaul
	if [[ "${_arg_menu[1]}" == "ptgaul" ]]; then

		_contigger_edges_gfa="${_arg_outdir}/result_3000/flye_cpONT/assembly_graph.gfa"
		_outdir="${_arg_outdir}/result_3000/flye_cpONT/ptdna"

		_mtcontigname="${_outdir}/mt.contig.name"
		_arg_unpolished_fasta="${_outdir}/circular_path_1_concatenated.fa"
		_arg_final_assembly="${_outdir}/pt.1.fa"

		if [[ "${_arg_menu[2]}" == "outfile" ]]; then
			_polap_log3_cmd mkdir -p "${_outdir}"
			gfatools stat -l "${_contigger_edges_gfa}" \
				>"${_outdir}/graph_final.txt" \
				2>"$_polap_output_dest"
			gfa_number_segments=$(grep "Number of segments:" "${_outdir}/graph_final.txt" | awk '{print $4}')

			if [[ "${gfa_number_segments}" -eq 3 ]]; then
				_polap_log0 "ptGAUL assembly's segments: ${gfa_number_segments}"
				_polap_log0 "you may use ptgaul menu:"
				_polap_log0 "  ptgaul 1"
				_polap_log0 "  ptgaul 1 1"
				_polap_log0 "  ptgaul 1 3"
				_polap_log0 "  ptgaul 2"
				_polap_log0 "  ptgaul 3"
				# mt.contig.name
				_polap_log0 "create file: ${_mtcontigname}"
				local string="edge_1,edge_2,edge_3"
				echo "$string" | tr ',' '\n' >"${_mtcontigname}"
				# extract ptDNA from ptGAUL's result
				_polap_log3_pipe "python \
          $script_dir/run-polap-py-find-plastid-gfa2fasta.py \
		        --gfa ${_contigger_edges_gfa} \
		        --seed ${_mtcontigname} \
		        --out ${_outdir} \
		        2>$_polap_output_dest"
				# polish one of the ptDNA sequences
				_polap_log0 "polishing ptDNA: ${_arg_unpolished_fasta}"
				_run_polap_polish
				_polap_log0 "ptGAUL polished assembly: ${_arg_final_assembly}"

			elif [[ "${gfa_number_segments}" -eq 1 ]]; then
				_polap_log0 "ptGAUL assembly's segments: ${gfa_number_segments}"
				_arg_unpolished_fasta="${_arg_outdir}/result_3000/ptGAUL_final_assembly/final_assembly.fasta"
				_polap_log0 "ptGAUL assembly: ${_arg_unpolished_fasta}"
				_polap_log0 "polishing ptDNA: ${_arg_unpolished_fasta}"
				_run_polap_polish
				_polap_log0 "ptGAUL polished assembly: ${_arg_final_assembly}"
			else
				_polap_log0 "ptGAUL assembly's segment counts: ${gfa_number_segments}"
				_polap_log0 "you may not use ptgaul menu at the moment."
				_polap_log0 "check: ${_contigger_edges_gfa}"
				_polap_log0 "create: ${_mtcontigname} with e.g., edge_1 that is originated from the ptDNA."
				_polap_log0 "then run disassemble ptgaul 2"
				_polap_log0 "then run disassemble ptgaul 3"
			fi
		fi

		# _mtcontigname
		if [[ "${_arg_menu[2]}" == "1" ]]; then
			if [[ "${_arg_menu[3]}" == "1" ]]; then
				_polap_log0 "create file: ${_mtcontigname}"
				local string="edge_1"
				echo "$string" | tr ',' '\n' >"${_mtcontigname}"
			elif [[ "${_arg_menu[3]}" == "3" ]]; then
				_polap_log0 "create file: ${_mtcontigname}"
				local string="edge_1,edge_2,edge_3"
				echo "$string" | tr ',' '\n' >"${_mtcontigname}"
			fi
		fi

		# extract
		if [[ "${_arg_menu[2]}" == "2" ]]; then
			_polap_log0 "extract ptDNA: ${_arg_unpolished_fasta}"
			_polap_log3_pipe "python \
          $script_dir/run-polap-py-find-plastid-gfa2fasta.py \
		        --gfa ${_contigger_edges_gfa} \
		        --seed ${_mtcontigname} \
		        --out ${_outdir} \
		        2>$_polap_output_dest"
		fi

		# polishing
		if [[ "${_arg_menu[2]}" == "3" ]]; then
			_polap_log0 "polishing ptDNA: ${_arg_unpolished_fasta}"
			_run_polap_polish
			_polap_log0 "ptGAUL polished assembly: ${_arg_final_assembly}"
		fi

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# menu: example
	if [[ "${_arg_menu[1]}" == "example" ]]; then
		if [[ "${_arg_menu[2]}" == "1" ]]; then
			"$0" x-ncbi-fetch-sra --sra SRR7153095
			"$0" x-ncbi-fetch-sra --sra SRR7161123
			"$0" get-mtdna --plastid --species "Eucalyptus pauciflora"
			cp o/00-bioproject/2-mtdna.fasta ref.fasta
		fi

		if [[ "${_arg_menu[2]}" == "2" ]]; then
			"$0" disassemble -l SRR7153095.fastq \
				-a SRR7161123_1.fastq -b SRR7161123_2.fastq
		fi

		if [[ "${_arg_menu[2]}" == "3" ]]; then
			"$0" disassemble -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq --disassemble-compare-to-fasta ref.fasta --disassemble-min-memory 9
		fi

		if [[ "${_arg_menu[2]}" == "4" ]]; then
			# stage 1: no short-read
			"$0" disassemble -o epauciflora -l SRR7153095.fastq --disassemble-s-max 100m --disassemble-compare-to-fasta ptdna-epauciflora.fa --disassemble-min-memory 9
		fi

		if [[ "${_arg_menu[2]}" == "5" ]]; then
			# stage 2: no short-read
			"$0" disassemble -o epauciflora -l SRR7153095.fastq -g 250k --disassemble-s 200m --disassemble-alpha 1.5 --disassemble-compare-to-fasta ptdna-epauciflora.fa --disassemble-min-memory 9 --contigger
		fi

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	# main
	#
	# select stages
	if [[ -z "${_arg_stages_include}" ]]; then
		_arg_stages_include="0-4"
		_arg_stages_exclude=""
	fi
	_include="${_arg_stages_include}"
	_exclude="${_arg_stages_exclude}" # Optional range or list of steps to exclude
	_stage_array=()
	_stage_array=($(_polap_parse_steps "${_include}" "${_exclude}"))

	_polap_log0 "Assemble plastid genomes by subsampling long-read data (${_arg_disassemble_i})"
	_polap_log1 "  input1: long-read: ${_arg_long_reads}"
	[[ -s "${_arg_long_reads}" ]] || return ${_POLAP_ERR_CMD_OPTION_LONGREAD}
	if [[ "${_arg_short_read1_is}" == "off" ]]; then
		_polap_log1 "  input2: no short-read1"
	else
		_polap_log1 "  input2: short-read1: ${_arg_short_read1}"
		[[ -s "${_arg_short_read1}" ]] || return ${_POLAP_ERR_CMD_OPTION_SHORTREAD}
	fi
	if [[ "${_arg_short_read2_is}" == "off" ]]; then
		_polap_log1 "  input3: no short-read2"
	else
		_polap_log1 "  input3: short-read2: ${_arg_short_read2}"
		[[ -s "${_arg_short_read2}" ]] || return ${_POLAP_ERR_CMD_OPTION_SHORTREAD}
	fi
	# Iteration over p x n for each i in --disassemble-i
	#
	# output: o/disassemble + --disassemble-i
	# output: o/disassemble/i
	_disassemble_i="${_disassemble_dir}/${_arg_disassemble_i}"
	_polap_log3_cmd mkdir -p "${_disassemble_i}"
	_polap_log1 "  output: ${_disassemble_i}"

	# output: o/disassemble/i/params.txt
	_disassemble_params_file="${_disassemble_i}/params.txt"
	(
		make_params_txt() {
			rm -f "${_disassemble_params_file}"
			printf "%s: %s\n" "I" "${_arg_disassemble_i}" >>"${_disassemble_params_file}"
			printf "%s: %s\n" "P" "${_arg_disassemble_p}" >>"${_disassemble_params_file}"
			printf "%s: %s\n" "N" "${_arg_disassemble_n}" >>"${_disassemble_params_file}"
			printf "%s: %s\n" "R" "${_arg_disassemble_r}" >>"${_disassemble_params_file}"
			printf "%s: %s\n" "A" "${_arg_disassemble_a}" >>"${_disassemble_params_file}"
		}
		make_params_txt
	)

	# Stage 0
	# short-read polishing preparation
	#
	# output: o/msbwt
	if _polap_contains_step 0 "${_stage_array[@]}"; then
		_run_polap_prepare-polishing
	fi

	# Stage 1
	# output: o/disassemble/i/1
	if _polap_contains_step 1 "${_stage_array[@]}"; then
		_disassemble_i_stage="${_disassemble_i}/1"
		_polap_log0 "  stage 1: assemble with a range of sample sizes"

		if [[ "${_arg_disassemble_c_is}" == "off" ]]; then
			_arg_contigger="on"
			_arg_polish="off"
		else
			_arg_contigger="off"
			_arg_polish="on"
		fi

		if [[ -z "${_arg_disassemble_s}" ]]; then
			_polap_log1 "  no option: --disassemble-s"
			_polap_log1 "    -> use a range of rates for subsampling long-read sequencing data"

			# short-read data
			if [[ "${_arg_short_read1_is}" == "on" ]] ||
				[[ "${_arg_short_read2_is}" == "on" ]]; then
				_polap_log1 "  use short-read data for genome size estimates"
				if [[ -z "${_arg_steps_include}" ]]; then
					if [[ "${_arg_disassemble_best}" == "off" ]]; then
						_arg_steps_include="1-16"
					else
						_arg_steps_include="11,12,13,14"
					fi
				fi
			else
				_polap_log1 "  no short-read data, so use a range of genome sizes"
				if [[ -z "${_arg_steps_include}" ]]; then
					_arg_steps_include="1-16"
					_arg_steps_exclude="2,3,6"
				fi
			fi
			_run_polap_step-disassemble 1
			_rstatus="$?"
			[[ "$_rstatus" -ne 0 ]] && return "$_rstatus"
		fi
	else
		# --disassemble-s -> only one sample or single s
		_polap_log1 "  --disassemble-s -> no iteration or no stage 1"
	fi

	# use summary1 to select one
	#
	if _polap_contains_step 2 "${_stage_array[@]}"; then
		_disassemble_i_stage="${_disassemble_i}/1"

		if [[ -z "${_arg_disassemble_s}" ]]; then
			_summary1_ordered="${_disassemble_i_stage}/summary1-ordered.txt"

			if [[ "${_arg_disassemble_best}" == "off" ]]; then
				_polap_log3_pipe "Rscript --vanilla $script_dir/run-polap-r-disassemble.R \
        --table ${_disassemble_i_stage}/summary1.txt \
        --out ${_disassemble_i_stage}/summary1-ordered.txt \
        --plot ${_disassemble_i_stage}/summary1-ordered.pdf \
        2>${_polap_output_dest}"
				# 2>&1"
				_rstatus="$?"
				if [[ "$_rstatus" -ne 0 ]]; then
					_polap_log0 "ERROR: run-polap-r-disassemble.R on ${_disassemble_i_stage}/summary1.txt"
				fi
			fi

			if [[ -s "${_summary1_ordered}" ]]; then
				# Use awk to extract the desired columns and save the output to a Bash variable
				output=$(awk -F'\t' 'NR==2 {print $1, $2, $4, $5, $9, $11}' "${_summary1_ordered}")
				n=$(awk 'END {print NR - 1}' "${_summary1_ordered}")
				if ((n < 5)); then
					_polap_log0 "  warning: the number of potential ptDNA assemblies is too small to select one"
					_polap_log0 "  suggestion: increase the subsample size --disassemble_p or the step size --disassemble_n"
					_polap_log0 "  suggestion: increase the subsample size --disassemble_b or the step size --disassemble_n"
				elif ((n < 2)); then
					_polap_log0 "ERROR: the number of potential ptDNA assemblies is less than 2"
					return "${_POLAP_ERR_SUBSAMPLE_TOO_FEW_CANDIDATES}"
				fi
				# Use read to split the output into individual variables
				read -r index size rate randomseed genomesize alpha <<<"$output"
				# Print the extracted variables
				_polap_log1 "  best Index: $index"
				_polap_log1 "  rate: $rate"
				_polap_log1 "  size: $size"
				_polap_log1 "  alpha: $alpha"
				_polap_log1 "  genome size: $genomesize"
				_arg_disassemble_s="$size"
				_arg_disassemble_alpha="$alpha"
				_arg_disassemble_n_is="on"
			else
				_polap_log0 "ERROR: no such file: ${_summary1_ordered}"
				_polap_log0 "ERROR: subsampling does not produce enough assemblies."
				die "SUGGESTION: increase the sampling size --disassemble-n"
			fi
		else
			_polap_log1 "  --disassemble-s -> no iteration or no stage 1"
		fi
	fi

	if [[ -n "${_arg_disassemble_stop_after}" ]]; then
		if [[ "${_arg_disassemble_stop_after}" == "stage1" ]]; then
			return 0
		fi
	fi

	_arg_steps_include=""
	if _polap_contains_step 3 "${_stage_array[@]}"; then
		_disassemble_i_stage="${_disassemble_i}/2"

		_arg_contigger="off"
		_arg_polish="on"

		# stage 2:
		_polap_log0 "  stage 2: assemble with a single set of parameters"
		if [[ "${_arg_disassemble_best}" == "off" ]]; then
			if [[ -z "${_arg_disassemble_s}" ]]; then
				_polap_log0 "ERROR: subsample size must be given either by the stage 1 or option --disassemble-s"
				die "ERROR: no --disassemble-s"
			fi
			if [[ -z "${_arg_disassemble_alpha}" ]]; then
				_polap_log0 "ERROR: subsampling minimum coverage alpha must be given either by the stage 1 or option --disassemble-s"
				die "ERROR: no --disassemble-alpha"
			fi
			# _arg_contigger="on"
			_arg_disassemble_a="${_arg_disassemble_s}"
			_arg_disassemble_b="${_arg_disassemble_s}"
			_arg_disassemble_b_is="on"
			_arg_disassemble_n="${_arg_disassemble_r}"
			_polap_log1 "  input1: --disassemble-s=${_arg_disassemble_s}"
		fi

		if [[ "${_arg_short_read1_is}" == "on" ]] ||
			[[ "${_arg_short_read2_is}" == "on" ]]; then
			_polap_log2 "  genome size estimated by the short-read data"
			if [[ -z "${_arg_steps_include}" ]]; then
				if [[ "${_arg_disassemble_best}" == "off" ]]; then
					_arg_steps_include="1-16"
					_arg_steps_exclude="9"
					check_file_existence "${_input_short_reads}"
				else
					_arg_steps_include="11,12,13,14"
				fi
			fi
		else
			if [[ -z "${_arg_steps_include}" ]]; then
				_arg_steps_include="1-16"
				_arg_steps_exclude="2,3,6,9"
			fi
			if [[ -z "${_arg_genomesize}" ]]; then
				_polap_log0 "ERROR: menu disassemble without short-read data and with --disassemble-s option"
				return "${_POLAP_ERR_CMD_OPTION_GENOMESIZE}"
			else
				_polap_log1 "  input2: genome size: ${_arg_genomesize}"
			fi
		fi

		_run_polap_step-disassemble 2
		_rstatus="$?"
		[[ "$_rstatus" -ne 0 ]] && return "$_rstatus"
	fi

	if _polap_contains_step 4 "${_stage_array[@]}"; then
		_disassemble_i_stage="${_disassemble_i}/2"
		_summary1_ordered="${_disassemble_i_stage}/summary1-ordered.txt"

		if [[ "${_arg_disassemble_best}" == "off" ]]; then
			_polap_log3_pipe "Rscript --vanilla $script_dir/run-polap-r-disassemble.R \
        --table ${_disassemble_i_stage}/summary1.txt \
        --out ${_disassemble_i_stage}/summary1-ordered.txt \
        --plot ${_disassemble_i_stage}/summary1-ordered.pdf \
        2>${_polap_output_dest}"
			# 2>&1"
		fi

		if [[ -s "${_summary1_ordered}" ]]; then
			# Use awk to extract the desired columns and save the output to a Bash variable
			output=$(awk -F'\t' 'NR==2 {print $1, $2, $4, $5, $9, $11}' "${_summary1_ordered}")
			n=$(awk 'END {print NR - 1}' "${_summary1_ordered}")
			if ((n < 5)); then
				_polap_log0 "  warning: the number of potential ptDNA assemblies is too small to select one"
				_polap_log0 "  suggestion: increase the replicate size --disassemble_r"
			elif ((n < 2)); then
				_polap_log0 "ERROR: the number of potential ptDNA assemblies is less than 2"
				return "${_POLAP_ERR_SUBSAMPLE_TOO_FEW_CANDIDATES}"
			fi
			# Use read to split the output into individual variables
			read -r index size rate randomseed genomesize alpha <<<"$output"
			# Print the extracted variables
			_polap_log1 "  best index: $index"
			_var_mtdna="${_disassemble_i_stage}/${index}/52-mtdna"
			_arg_unpolished_fasta="${_var_mtdna}/ptdna.0.fa"
			cp -p "${_var_mtdna}/ptdna.1.fa" "${_disassemble_i}/pt.1.fa"
		else
			die "ERROR: no such file: ${_summary1_ordered}"
		fi
	fi

	# Disable debugging if previously enabled
	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	[ "$DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_polap_step-disassemble {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	_arg_plastid="on"
	local i
	local j
	local col1 col2 col3
	local _command1
	local _disassemble_dir
	local _input_short_reads
	local _include
	local _exclude
	local _datasize_array
	local _datasize_array_length
	local _i_last
	local _index_table
	local sampling_datasize
	local genomesize
	local _size_index
	local _alpha
	local long_total
	local short_total
	local long_rate_sample
	local short_rate_sample
	local long_sample_seed
	local _pre_alpha
	local _outdir
	local _sampling_datasize_bp
	local _total_length
	local _total_short_read
	local _short_read1
	local _long_read
	local is_stop
	local _summary1_file
	local resume=""
	local num_circular_paths
	local num_circular_nodes
	local _expected_genome_size
	local _contigger_edges_gfa
	local _contigger_edges_fasta
	local peak_ram_size_gb
	local peak_ram
	local _outdir_draft_assembly_size
	local _draft_assembly_size
	local _draft_assembly_size_bp
	local _sampling_datasize_bp
	local gfa_number_segments
	local gfa_number_links
	local gfa_total_segment_length
	local _ga_annotation_all
	local _mtcontigname
	local _ptdna
	local _ptdir
	local ptdna_file
	local nseq
	local pident
	local gc
	local length
	local coverage
	local _total_iterations
	local time_per_iteration
	local remaining_iterations
	local remaining_time
	local _terminal_width
	local status
	local _var_mtdna
	local _restarted_fasta
	local _circular_path_fasta
	local _folder_with_max_coverage
	local _step_array

	source "$script_dir/polap-variables-common.sh"
	source "$script_dir/polap-function-disassemble-seeds.sh"

	help_message=$(
		cat <<HEREDOC
Plastid genome assembly step-by-step using long-read data without reference

Inputs
------

- long-read data: ${_arg_long_reads}
- short-read data 1: ${_arg_short_read1}
- short-read data 2 (optional): ${_arg_short_read2}

Outputs
-------

- plastid genome assembly
- trace plots for the features of plastid genome assemblies
- ${_arg_outdir}/disassemble/<number>

Arguments
---------

Stage 1:
-o ${_arg_outdir}
-l ${_arg_long_reads}: a long-read fastq data file
-a ${_arg_short_read1}: a short-read fastq data file 1
-b ${_arg_short_read2} (optional): a short-read fastq data file 2
--disassemble-min-memory ${_arg_disassemble_min_memory}: the minimum memory in Gb
--disassemble-a ${_arg_disassemble_a}: the smallest long read step size
--disassemble-b ${_arg_disassemble_b}: the largest long read step size
--disassemble-n ${_arg_disassemble_n}: the number of steps
--disassemble-alpha ${_arg_disassemble_alpha}: the minimum Flye's disjointig coverage
--disassemble-delta ${_arg_disassemble_delta}: the move size of alpha
--disassemble-compare-to-fasta <FASTA>
--start-index <INDEX>: from <INDEX>
--end-index <INDEX>: to <INDEX> - 1
--steps-include <STPES>: STEPS can be 1,2,3 or 1-15
--steps-exclude <STPES>: STEPS can be 1,2,3 or 1-15
-t ${_arg_threads}: the number of CPU cores
--random-seed <arg>: 5-digit number
--species "species name"
--no-contigger for a single final step: use the long-read polished gfa

--disassemble-m ${_arg_disassemble_m}: the threshold to change alpha
--disassemble-s-max: the threshold to stop the iteration per genome size

Stage 2:
--disassemble-s: the threshold to change alpha
--disassemble-alpha ${_arg_disassemble_alpha}: the minimum Flye's disjointig coverage

Menus
-----

- help: display this help message
- redo: overwrite previous results
- view: show some results
- reset or clean: delete output files 

Usages
------
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads} -a ${_arg_short_read1}
$(basename "$0") ${_arg_menu[0]} -l ${_arg_long_reads} -a ${_arg_short_read1} -b ${_arg_short_read2}
$0 example 1
$0 example 2
tail -f o/polap.log
$0 example 3

Examples
--------
1:
$(basename "$0") x-ncbi-fetch-sra --sra SRR7153095
$(basename "$0") x-ncbi-fetch-sra --sra SRR7161123
$(basename "$0") get-mtdna --plastid --species "Eucalyptus pauciflora"
cp o/00-bioproject/2-mtdna.fasta ref.fasta
2:
$(basename "$0") disassemble -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq --disassemble-compare-to-fasta ref.fasta --disassemble-min-memory 9
$(basename "$0") disassemble view
3:
$(basename "$0") disassemble view 2
$(basename "$0") disassemble -o o2 --anotherdir o --disassemble-s 314588465 --disassemble-alpha 3.5 -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq

For ptDNA of Juncus validus
$(basename "$0") get-mtdna -o jvalidus --plastid --species "Juncus validus"
$(basename "$0") disassemble -o jvalidus -l SRR21976089.fastq -a SRR21976091_1.fastq -b SRR21976091_2.fastq --disassemble-min-memory 9 -v --disassemble-a 50mb --disassemble-n 20
src/polap.sh disassemble -l SRR21976089.fastq -a SRR21976091_1.fastq -b SRR21976091_2.fastq -o jvalidus --disassemble-min-memory 9 -v --disassemble-a 50mb --disassemble-n 20
HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return 0
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		case "${_arg_menu[2]}" in
		1)
			_polap_log0_column "${_disassemble_dir}/summary1.txt"
			;;
		2)
			_polap_log0_column "${_disassemble_dir}/summary2.txt"
			;;
		*)
			_polap_log0_column "${_disassemble_dir}/summary1.txt"
			;;
		esac

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "reset" ]] ||
		[[ "${_arg_menu[1]}" == "clean" ]]; then
		case "${_arg_menu[2]}" in
		1)
			rm -f "${_polap_var_outdir_long_total_length}"
			;;
		2)
			rm -f "${_input_short_reads}"
			;;
		3)
			rm -f "${_polap_var_outdir_short_total_length}"
			;;
		4)
			rm -rf "${_disassemble_dir}"
			;;
		*)
			echo default
			;;
		esac

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	if [[ "${_arg_menu[1]}" == "example" ]]; then
		if [[ "${_arg_menu[2]}" == "1" ]]; then
			"$0" x-ncbi-fetch-sra --sra SRR7153095
			"$0" x-ncbi-fetch-sra --sra SRR7161123
			"$0" get-mtdna --plastid --species "Eucalyptus pauciflora"
			cp o/00-bioproject/2-mtdna.fasta ref.fasta
		fi
		if [[ "${_arg_menu[2]}" == "2" ]]; then
			"$0" disassemble -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq --disassemble-compare-to-fasta ref.fasta --disassemble-min-memory 9
			"$0" disassemble view
		fi
		if [[ "${_arg_menu[2]}" == "3" ]]; then
			"$0" disassemble view 2
			"$0" disassemble -o o2 --anotherdir o --disassemble-s 314mb --disassemble-alpha 3.5 -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq
		fi

		# Disable debugging if previously enabled
		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		[ "$DEBUG" -eq 1 ] && set +x
		return 0
	fi

	i=0
	gfa_number_segments=0
	gfa_number_links=0
	gfa_total_segment_length=0
	num_circular_paths=0
	num_circular_nodes=0
	length=0
	gc=0
	coverage=0
	_input_short_reads="${_arg_outdir}/s1.fq"
	_disassemble_dir="${_arg_outdir}/disassemble/${_arg_disassemble_i}/${1}"
	_index_table="${_disassemble_dir}/index.txt"
	_short_read1=${_disassemble_dir}/s.fq
	_long_read=${_disassemble_dir}/l.fq.gz
	_summary1_file="${_disassemble_dir}/summary1.txt"

	_polap_log1 "  assemble a plastid genome without reference"
	if [[ "${1}" -eq 1 ]] || [[ "${1}" -eq 2 ]]; then
		_polap_log2 "    disassemble stage ${1}."
	else
		_polap_log0 "ERROR: disassemble stage must be 1 or 2."
	fi
	_polap_log2 "    input1: long-read: ${_arg_long_reads}"
	if [[ "${_arg_short_read1_is}" == "off" ]]; then
		_polap_log2 "    input2: no short-read1"
	else
		_polap_log2 "    input2: short-read1: ${_arg_short_read1}"
	fi
	if [[ "${_arg_short_read2_is}" == "off" ]]; then
		_polap_log2 "    input3: no short-read2"
	else
		_polap_log2 "    input3: short-read2: ${_arg_short_read2}"
	fi
	_polap_log2 "    input4: option --plastid: ${_arg_plastid}"

	if [[ -n "${_arg_species}" ]]; then
		_polap_log2 "  input6: species: ${_arg_species}"
		_run_polap_get-mtdna
		_arg_disassemble_c="${_arg_outdir}/00-bioproject/2-mtdna.fasta"
	elif [[ -n "${_arg_disassemble_c}" ]]; then
		_polap_log2 "  input6: compare to fasta: ${_arg_disassemble_c}"
		check_file_existence "${_arg_disassemble_c}"
	else
		_polap_log2 "  input6: neither species nor fasta file is given"
		_polap_log2 "    inference of ptDNA"
		_polap_log2 "    one of the sequences of an assembly graph is chosen"
	fi
	_polap_log2 "  output: ${_disassemble_dir}"

	# select steps
	local _include="${_arg_steps_include}"
	local _exclude="${_arg_steps_exclude}" # Optional range or list of steps to exclude
	_step_array=()

	_step_array=($(_polap_parse_steps "${_include}" "${_exclude}"))
	_rstatus="$?"
	if [[ "${_rstatus}" -ne 0 ]]; then
		_polap_log0 "ERROR: Error parsing steps."
		return "${_POLAP_ERR_CMD_OPTION_STEPS}"
	fi

	_polap_log1 "  executing steps: ${_step_array[*]}"

	# steps
	# 1. long-read count
	# 2. short-read concatenate
	# 3. short-read count
	# 4. loop index
	# 5. long-read sample
	# 6. short-read sample
	# 7. genome size estimate
	# 8. flye run
	# 9. alpha adjustment
	# 10. assembly summary
	# 11. annotate
	# 12. seeds
	# 13. ptDNA extraction
	# 14. best ptDNA
	# 15. ptDNA evaluation
	# 16.

	# 1. long-read count
	if _polap_contains_step 1 "${_step_array[@]}"; then
		_polap_log1 "  step 1: count the total number of bases in the long-read data"
		_polap_log2 "    input1: long-read: ${_long_read}"
		(
			step1() {
				if [ -s "${_polap_var_outdir_long_total_length}" ] &&
					[ "${_arg_redo}" = "off" ]; then
					_polap_log2 "  found: ${_polap_var_outdir_long_total_length}, skipping total-length-long ..."
				else
					_run_polap_total-length-long
				fi
			}
			step1
		)
	fi

	# 2. short-read concatenate
	# combine short-read dataset files
	# create o/s1.fq to deal with a single or two short-read files
	# overwrite the s1.fq
	if _polap_contains_step 2 "${_step_array[@]}"; then
		_polap_log1 "  step 2: concatenate the short-read data"
		_polap_log2 "    output: ${_input_short_reads}"
		if { [[ "${_arg_short_read1_is}" == "on" ]] ||
			[[ "${_arg_short_read2_is}" == "on" ]]; } &&
			[[ "${_arg_disassemble_best}" == "off" ]]; then
			(
				step2() {
					if [ -s "${_input_short_reads}" ] && [ "${_arg_redo}" = "off" ]; then
						_polap_log2 "  found: ${_input_short_reads}, skipping ..."
					else
						if [[ "${_arg_short_read1_is}" == "on" ]] ||
							[[ "${_arg_short_read2_is}" == "on" ]]; then
							_polap_concatenate_fastq_to "${_input_short_reads}" "${_arg_short_read1}" "${_arg_short_read2}"
						elif [[ "${_arg_short_read2_is}" == "on" ]]; then
							_polap_concatenate_fastq_to "${_input_short_reads}" "${_arg_short_read2}"
						else
							_polap_concatenate_fastq_to "${_input_short_reads}" "${_arg_short_read1}"
						fi
					fi
				}
				step2
			)
		fi
		check_file_existence "${_input_short_reads}"
	fi

	# 3. short-read count
	if _polap_contains_step 3 "${_step_array[@]}"; then
		_polap_log1 "  step 3: count the total number of bases in the short-read data"
		_polap_log2 "    input1: ${_input_short_reads}"
		check_file_existence "${_input_short_reads}"
		if { [[ "${_arg_short_read1_is}" == "on" ]] ||
			[[ "${_arg_short_read2_is}" == "on" ]]; } &&
			[[ "${_arg_disassemble_best}" == "off" ]]; then
			(
				step3() {
					if [ -s "${_polap_var_outdir_short_total_length}" ] && [ "${_arg_redo}" = "off" ]; then
						_polap_log2 "  found: ${_polap_var_outdir_short_total_length}, skipping total-length-short ..."
					else
						_polap_total-length-short \
							"${_input_short_reads}" \
							"${_polap_var_outdir_short_total_length}"

						local _total_length=$(<"${_polap_var_outdir_short_total_length}")
						local _total_short_read=$(_polap_utility_convert_bp ${_total_length})
						_polap_log2 "    the size of the short-read dataset: ${_total_short_read}"
					fi
				}
				step3
			)
		fi
	fi

	# 4. loop index
	if _polap_contains_step 4 "${_step_array[@]}"; then
		_polap_log1 "  step 4: prepare for the loop iteration"
		if [ -d "${_disassemble_dir}" ] && [ "${_arg_redo}" = "on" ]; then
			_polap_log3_cmd rm -rf "${_disassemble_dir}"
		fi
		_polap_log3_cmd mkdir -p "${_disassemble_dir}"

		if [[ "${_arg_disassemble_b_is}" == "on" ]]; then
			_run_polap_step-disassemble-make-index-for-b "${_index_table}"
		else
			_run_polap_step-disassemble-make-index-for-p "${_index_table}"
		fi
		_alpha="${_arg_disassemble_alpha}"

		# Iterate using indices of the data size
		_datasize_array=()
		while IFS=$'\t' read -r col1 col2 col3; do
			_datasize_array+=("$col1 $col2 $col3") # Append each row as a string
		done <"${_disassemble_dir}/index.txt"
		_datasize_array_length=${#_datasize_array[@]}
		_total_iterations="${_datasize_array_length}"

		_size_index="${#_datasize_array[@]}"
		if ((_arg_end_index > 0)); then
			_size_index="${_arg_end_index}"
		fi

		# Summary
		#
		# long_total
		# long_rate_sample
		# long_sample_seed
		#
		# short_total
		# short_rate_sample
		# short_sample_seed
		#
		# _expected_genome_size
		#
		# peak_ram_size_gb
		#
		# alpha
		#
		# gfa_number_segments
		# gfa_number_links
		# gfa_total_segment_length
		#
		# note: annotation result e.g., number of MT/PT genes
		# note: seeds e.g., number of seeds
		#
		# num_circular_paths
		# num_circular_nodes
		#
		# coverage
		#
		# pident
		#
		# length
		# gc

		# Add header to the output file
		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
			"index" \
			"sampling_datasize" \
			"long_total" "long_rate_sample" "long_sample_seed" \
			"short_total" "short_rate_sample" "short_sample_seed" \
			"expected_genome_size" \
			"peak_ram_size_gb" \
			"alpha" \
			"gfa_number_segments" "gfa_number_links" "gfa_total_segment_length" \
			"num_circular_paths" "num_circular_nodes" \
			"time" \
			"coverage_ref" \
			"coverage_target" \
			"pident" \
			"length" \
			"gc" \
			"j_candidate" \
			>"${_summary1_file}"
		is_stop="off"
	else
		rm -f "${_index_table}"
		>"${_index_table}"
		_alpha="${_arg_disassemble_alpha}"
		_total_iterations=0
		is_stop="on"
	fi

	_i_last=0
	_terminal_width=$(tput cols) # Get terminal width
	cat <"${_index_table}" | while IFS=$'\t' read -r i sampling_datasize genomesize; do
		if [[ "${is_stop}" == "on" ]]; then
			break
		fi

		if ((i < _arg_start_index)); then
			continue
		fi

		if ((i >= _size_index)); then
			break
		fi

		if [[ -n "${_arg_disassemble_s_max}" ]]; then
			if [[ "${sampling_datasize}" -gt "${_arg_disassemble_s_max}" ]]; then
				break
			fi
		fi

		_polap_set_start_time
		_i_last=$((_i_last + 1))
		_pre_alpha="${_alpha}"
		_outdir=${_disassemble_dir}/$i
		_contigger_edges_gfa="${_outdir}/30-contigger/graph_final.gfa"
		_contigger_edges_fasta="${_outdir}/30-contigger/graph_final.fasta"
		_mtcontigname="${_outdir}/mt.contig.name"
		_ga_annotation_all="${_outdir}/assembly_info_organelle_annotation_count-all.txt"
		_var_mtdna="${_outdir}/52-mtdna"
		_ptdna="${_var_mtdna}/circular_path_1_concatenated.fa"

		_sampling_datasize_bp=$(_polap_utility_convert_bp "${sampling_datasize}")
		_genomesize_bp=$(_polap_utility_convert_bp "${genomesize}")
		_polap_log2 "  (i=${i}):S=${_sampling_datasize_bp}:G=${_genomesize_bp}:alpha=${_pre_alpha}"
		_polap_log3_cmd mkdir -p "${_outdir}"

		# 5. Sample long-read data
		#
		long_total=-1
		long_rate_sample=-1
		long_sample_seed=-1
		if _polap_contains_step 5 "${_step_array[@]}"; then
			_polap_log1 "  step 5: subsample long-read data"
			long_total=$(<"${_polap_var_outdir_long_total_length}")
			if [[ -s "${_long_read}" ]] && [[ "${_arg_test}" == "on" ]]; then
				_polap_log2 "  (test only) found: ${_long_read}, and skipping ..."
			else
				long_rate_sample=$(echo "scale=10; ${sampling_datasize}/${long_total}" | bc)
				long_sample_seed=${_arg_random_seed:-$RANDOM}
				_polap_log2 "    input1: sampling rate: ${long_rate_sample}"
				_polap_log2 "    input2: random seed: ${long_sample_seed}"
				_polap_log2 "    input3: long read data: ${_arg_long_reads}"
				_polap_log2 "    random seed: ${long_sample_seed}"
				_polap_log2 "    output: ${_long_read}"
				_polap_log3_pipe "seqkit sample \
          -p ${long_rate_sample} \
		      -s ${long_sample_seed} \
		      ${_arg_long_reads} \
		      -o ${_long_read} 2>${_polap_output_dest}"
			fi
		fi

		# 6. Sample short-read data
		# We utilize the identical sampling rate employed
		# by the long-read sampling method, which is set to default.
		#
		short_total=-1
		short_rate_sample=-1
		short_sample_seed=-1
		if _polap_contains_step 6 "${_step_array[@]}"; then
			check_file_existence "${_input_short_reads}"
			_polap_log1 "  step 6: subsample the short-read data"
			short_total=$(<"${_polap_var_outdir_short_total_length}")
			if [[ -s "${_short_read1}" ]] && [[ "${_arg_test}" == "on" ]]; then
				_polap_log2 "  (test only) found: ${_short_read1}, and skipping ..."
			else
				if [[ "${_arg_disassemble_b_short}" == "on" ]]; then
					_polap_log1 "    use sample size: ${sampling_datasize} bp"
					short_rate_sample=$(echo "scale=10; ${sampling_datasize}/${short_total}" | bc)
				else
					_polap_log1 "    use sample rate: ${long_rate_sample}"
					short_rate_sample="${long_rate_sample}"
				fi
				short_sample_seed=${_arg_random_seed:-$RANDOM}
				_polap_log2 "    input1: sampling rate: ${short_rate_sample}"
				_polap_log2 "    input2: random seed: ${short_sample_seed}"
				_polap_log2 "    input3: short read data: ${_input_short_reads}"
				_polap_log2 "    random seed: ${short_sample_seed}"
				_polap_log2 "    output: ${_short_read1}"
				_polap_log3_pipe "seqkit sample \
          -p ${short_rate_sample} \
          -s ${short_sample_seed} \
          ${_input_short_reads} \
          -o ${_short_read1} 2>${_polap_output_dest}"
			fi
		fi

		# 7. Estitame the genome assembly size
		#
		_expected_genome_size=-1
		if _polap_contains_step 7 "${_step_array[@]}"; then
			_polap_log1 "  step 7: genome size using the subsample of the short-read data"
			if [[ "${_arg_short_read1_is}" == "off" ]] &&
				[[ "${_arg_short_read2_is}" == "off" ]]; then
				# Check if the genomesize is a positive number and not empty
				if [[ -n $genomesize && $genomesize =~ ^[0-9]+$ ]]; then
					_expected_genome_size="${genomesize}"
					_polap_log2 "    output: given genome size: ${_expected_genome_size} bp"
				else
					die "$genomesize is NOT a valid positive number or is empty."
				fi
			else
				if [[ -s "${_outdir}/short_expected_genome_size.txt" ]] &&
					[[ "${_arg_test}" == "on" ]]; then
					_polap_log2 "  found: ${_outdir}/short_expected_genome_size.txt, and skipping ..."
				else
					check_file_existence "${_short_read1}"
					_polap_log2 "    input1: ${_short_read1}"
					_polap_log2 "    output: ${_outdir}/short_expected_genome_size.txt"
					_polap_find-genome-size \
						"${_short_read1}" \
						"${_outdir}"
				fi
				_expected_genome_size=$(<"${_outdir}/short_expected_genome_size.txt")
				_polap_log2 "    output: estimated genome size: ${_expected_genome_size} bp"
			fi
		fi

		# 8. flye
		#
		peak_ram_size_gb=-1
		peak_ram=-1
		if _polap_contains_step 8 "${_step_array[@]}"; then
			_polap_log1 "  step 8: cFlye on the subsample of the long-read data"
			# no resume feature
			if [[ -d "${_outdir}/00-assembly" ]]; then
				_polap_log2 "  found: ${_outdir}/00-assembly, and skipping ..."
			else
				# cflye assemble stage
				_polap_log1 "  executing whole-genome assembly using cflye ... be patient!"
				_polap_log2 "    input1: long-read: ${_long_read}"
				_polap_log2 "    input2: genome size: ${_expected_genome_size}"
				_polap_log2 "    input3: alpha: ${_alpha}"
				_polap_log2 "    outdir: ${_outdir}"
				_run_polap_step-disassemble-cflye \
					"${_outdir}" "${_long_read}" "${_expected_genome_size}" "${_alpha}"
			fi

			# check the peak memory of cflye
			if [[ -d "${_outdir}/30-contigger" ]] &&
				[[ -s "${_contigger_edges_fasta}" ]]; then
				peak_ram_size_gb=$(grep "Peak RAM usage" "${_outdir}/cflye.log" | awk 'NR==2 {print $7}')
				peak_ram=$(grep "Peak RAM usage" "${_outdir}/cflye.log" | awk 'NR==2 {print $7, $8}')
				if [[ -z "${peak_ram}" ]]; then
					die "ERROR: cflye failed to run."
				fi
				if ((peak_ram_size_gb < _arg_disassemble_min_memory)); then
					_polap_log2 "  cFlye: used memory is less than ${_arg_disassemble_min_memory} Gb ($i): $peak_ram"
				else
					_polap_log2 "  cFlye: used memory is not less than ${_arg_disassemble_min_memory} Gb ($i): $peak_ram"
					_polap_log2 "  exit the disassemble menu."
					is_stop="on"
				fi

				# copy assembly_graph.gfa to graph_final.gfa if --no-contigger
				#
				if [[ "${_arg_contigger}" == "off" ]]; then
					if [[ -s "${_outdir}/assembly_graph.gfa" ]]; then
						if [[ ! -s "${_contigger_edges_gfa}.backup" ]]; then
							_polap_log2 "  backup ${_contigger_edges_gfa}"
							_polap_log3_cmd cp -p "${_contigger_edges_gfa}" "${_contigger_edges_gfa}.backup"
						fi
						_polap_log2 "  use ${_outdir}/assembly_graph.gfa by copying graph_final.gfa"
						_polap_log3_cmd cp -p "${_outdir}/assembly_graph.gfa" "${_contigger_edges_gfa}"
					else
						_polap_log0 "ERROR: cFlye did not finish with the polishing step."
						die "ERROR: no such file: ${_outdir}/assembly_graph.gfa"
					fi
				fi
			else
				_polap_log2 "  cFlye: no genome assembly"
			fi
		fi

		# 9. Adjust the alpha
		#
		_outdir_draft_assembly_size="${_outdir}/draft_assembly_size.txt"
		if _polap_contains_step 9 "${_step_array[@]}"; then
			_polap_log1 "  step 9: adjust the alpha"
			_polap_log2 "    input1: ${_outdir}/00-assembly/draft_assembly.fasta"
			_polap_log2 "    output: ${_outdir_draft_assembly_size}"

			# get the total number of bases in the flye draft assembly sequences
			_polap_log3_cmd rm -rf "${_outdir_draft_assembly_size}"
			if [[ -s "${_outdir}/00-assembly/draft_assembly.fasta" ]]; then
				_polap_log3_pipe "seqkit stats -Ta ${_outdir}/00-assembly/draft_assembly.fasta |\
		      csvtk cut -t -f sum_len |\
		      csvtk del-header \
		      >${_outdir_draft_assembly_size}"
				_draft_assembly_size=$(<"${_outdir_draft_assembly_size}")
			else
				_draft_assembly_size=0
			fi

			# adjust alpha
			if ((_arg_disassemble_m < _draft_assembly_size)); then
				_alpha=$(_polap_adjust_alpha_by_delta "${_alpha}" \
					"${_arg_disassemble_delta}" \
					"increase")
				_polap_log1 "    alpha increased: ${_alpha}"
			else
				_alpha=$(_polap_adjust_alpha_by_delta "${_alpha}" \
					"${_arg_disassemble_delta}" \
					"decrease")
				_polap_log1 "    alpha decrease (if not negative): ${_alpha}"
			fi
		fi

		# 10. Summary of assembly
		#
		gfa_number_segments=-1
		gfa_number_links=-1
		gfa_total_segment_length=-1
		if _polap_contains_step 10 "${_step_array[@]}"; then
			_polap_log1 "  step 10: summary statistics of gfa assembly graph"
			if [[ -s "${_contigger_edges_gfa}" ]] &&
				[[ -s "${_contigger_edges_fasta}" ]]; then
				_polap_log2 "    input1: ${_contigger_edges_gfa}"
				_polap_log2 "    output: ${_outdir}/graph_final.txt"
				_polap_log3_pipe "gfatools stat -l ${_contigger_edges_gfa} \
          >${_outdir}/graph_final.txt 2>$_polap_output_dest"
				gfa_number_segments=$(grep "Number of segments:" "${_outdir}/graph_final.txt" | awk '{print $4}')
				gfa_number_links=$(grep "Number of links:" "${_outdir}/graph_final.txt" | awk '{print $4}')
				gfa_total_segment_length=$(grep "Total segment length:" "${_outdir}/graph_final.txt" | awk '{print $4}')
			else
				if [[ ! -s "${_contigger_edges_gfa}" ]]; then
					_polap_log2 "    no input: ${_contigger_edges_gfa}"
				fi
				if [[ ! -s "${_contigger_edges_fasta}" ]]; then
					_polap_log2 "    no input: ${_contigger_edges_fasta}"
				fi
			fi
		fi

		# 11. Annotate
		#
		if _polap_contains_step 11 "${_step_array[@]}"; then
			_polap_log1 "  step 11: annotate to find the contig with more plastid genes than mitochondrial"
			if [[ -s "${_contigger_edges_gfa}" ]] &&
				[[ -s "${_contigger_edges_fasta}" ]]; then
				_polap_log2 "    input1: ${_contigger_edges_gfa}"
				_polap_log2 "    output: ${_ga_annotation_all}"
				if [[ -s "${_ga_annotation_all}" ]]; then
					_polap_log2 "  found: ${_ga_annotation_all}, so skipping mt/pt annotation ..."
				else
					polap_annotate "${_contigger_edges_gfa}" "${_ga_annotation_all}"
				fi
			else
				if [[ ! -s "${_contigger_edges_gfa}" ]]; then
					_polap_log2 "    no input: ${_contigger_edges_gfa}"
				fi
				if [[ ! -s "${_contigger_edges_fasta}" ]]; then
					_polap_log2 "    no input: ${_contigger_edges_fasta}"
				fi
			fi
		fi

		# 12. seeds
		# create mt.contig.name.txt
		# TODO:
		# filter gfa using depth from _mtcontigname.
		# filter out contigs with depths < smallest of mtcontig * 0.2
		if _polap_contains_step 12 "${_step_array[@]}"; then
			_polap_log1 "  step 12: seeds"
			if [[ -s "${_contigger_edges_gfa}" ]] &&
				[[ -s "${_contigger_edges_fasta}" ]] &&
				[[ -s "${_ga_annotation_all}" ]]; then
				_polap_log2 "    input1: ${_contigger_edges_gfa}"
				_polap_log2 "    output: ${_mtcontigname}"
				# depth-filtering
				if [[ -s "${_mtcontigname}" ]]; then
					_polap_log2 "  found: ${_mtcontigname}, so skipping seeds ..."
				else
					polap_disassemble-seeds "${_contigger_edges_gfa}" \
						"${_ga_annotation_all}" \
						"${_mtcontigname}"
				fi
			else
				if [[ ! -s "${_contigger_edges_gfa}" ]]; then
					_polap_log2 "    no input: ${_contigger_edges_gfa}"
				fi
				if [[ ! -s "${_contigger_edges_fasta}" ]]; then
					_polap_log2 "    no input: ${_contigger_edges_fasta}"
				fi
				if [[ ! -s "${_ga_annotation_all}" ]]; then
					_polap_log2 "    no input: ${_ga_annotation_all}"
				fi
			fi
		fi

		# 13. Extract ptDNA sequences
		#
		# locate the circular paths
		# extract candidate circular ptDNA sequences
		num_circular_paths=-1
		num_circular_nodes=-1
		if _polap_contains_step 13 "${_step_array[@]}"; then
			_polap_log1 "  step 13: extract ptDNA sequences from the gfa assembly graph"
			if [[ -s "${_contigger_edges_gfa}" ]] &&
				[[ -s "${_contigger_edges_fasta}" ]] &&
				[[ -s "${_mtcontigname}" ]]; then
				_polap_log2 "    input1: ${_contigger_edges_gfa}"
				_polap_log2 "    input2: ${_mtcontigname}"
				_polap_log2 "    outdir: ${_var_mtdna}"
				_polap_log3_cmd rm -rf "${_var_mtdna}"
				_polap_log3_cmd mkdir -p "${_var_mtdna}"
				_polap_log3_pipe "python \
          $script_dir/run-polap-py-find-plastid-gfa2fasta.py \
		        --gfa ${_contigger_edges_gfa} \
		        --seed ${_mtcontigname} \
		        --out ${_var_mtdna} \
		        2>$_polap_output_dest"
				if [[ -s "${_var_mtdna}/circular_path_count.txt" ]]; then
					num_circular_paths=$(<"${_var_mtdna}/circular_path_count.txt")
					num_circular_nodes=$(<"${_var_mtdna}/circular_path_nodes.txt")
					_polap_log3_cat "${_var_mtdna}/circular_path.txt"
					if [[ "${num_circular_paths}" -eq 4 ]] ||
						[[ "${num_circular_paths}" -eq 2 ]]; then
						_polap_log2 "  circular_path count is 2 or 4."
					else
						_polap_log2 "  circular_path count is neither 2 nor 4."
						_polap_log3_cmd rm -f "${_var_mtdna}"/circular_path_*_concatenated.fa
					fi
				fi
			else
				if [[ -s "${_contigger_edges_gfa}" ]]; then
					_polap_log2 "    no input: ${_contigger_edges_gfa}"
				fi
				if [[ -s "${_contigger_edges_fasta}" ]]; then
					_polap_log2 "    no input: ${_contigger_edges_fasta}"
				fi
				if [[ -s "${_mtcontigname}" ]]; then
					_polap_log2 "    no input: ${_mtcontigname}"
				fi
			fi
		fi

		# 14. Pick one sequence based on the reference, if any.
		# Use the first sequence if there is no such reference.
		#
		# for each ptDNA of the potential ptDNA sequences,
		#   select one candidate ptDNA
		#   rearrange it so that we could do a pairwise sequence alignment
		_coverage_ref=-1
		_coverage_target=-1
		_j_candidate=-1
		_circular_path_fasta="${_var_mtdna}/circular_path_1_concatenated.fa"
		_arg_unpolished_fasta="${_var_mtdna}/ptdna.0.fa"
		_arg_final_assembly="${_var_mtdna}/ptdna.1.fa"
		if _polap_contains_step 14 "${_step_array[@]}"; then
			_polap_log1 "  step 14: choose one of the multiple candidate ptDNA sequences"

			if { [[ -s "${_arg_unpolished_fasta}" ]] ||
				[[ -s "${_arg_final_assembly}" ]]; } &&
				[[ "${_arg_test}" == "on" ]]; then
				_polap_log2 "  found: ${_arg_unpolished_fasta}, so skipping selecting one ptDNA ..."
				_polap_log2 "  found: ${_arg_final_assembly}, so skipping selecting one ptDNA ..."
			fi

			num_circular_paths_file_count=$(find "${_var_mtdna}" -maxdepth 1 -name "circular_path_*_concatenated.fa" 2>/dev/null | wc -l)
			if [[ "${num_circular_paths_file_count}" -eq 4 ]] ||
				[[ "${num_circular_paths_file_count}" -eq 2 ]]; then
				_polap_log2 "    circular_path count is 2 or 4."
				_polap_log2 "    input1: ${_circular_path_fasta}"
				_polap_log2 "    output1: unpolished: ${_arg_unpolished_fasta}"
				_polap_log2 "    output2: polished: ${_arg_final_assembly}"

				if [[ -n "${_arg_disassemble_c}" ]]; then
					_polap_log1 "  step 14-1: use alignment coverage to select the best ptDNA using reference: ${_arg_disassemble_c}"
					_polap_log2 "    input1: ${_arg_disassemble_c}"
					_polap_log2 "    input2: e.g., ${_circular_path_fasta}"
					_polap_log2 "    outdir: ${_var_mtdna}"

					# 14-1. compute the coverage for each candidate
					# we typicially could have 4 or 2 candidates.
					# we might have more than 4 candidates, which might not be a good
					# situation.
					# FIXME: coverage needs some fix.
					# _coverage_ref: how much parts of the reference is coveraged in the alignment?
					# _coverage_target: how much parts of the target sequence?
					# both of these two coverage values need to be near 1.
					# we used only the _coverage_ref now.
					shopt -s nullglob # Enable nullglob so that the loop gets no elements if no match
					for _ptdna in "${_var_mtdna}"/circular_path_*_concatenated.fa; do
						j=$(_polap_utility_extract_number_circular_path "${_ptdna}")
						_ptdir="${_var_mtdna}/${j}"
						_polap_log3_pipe "python $script_dir/run-polap-py-compare2ptdna.py \
    		      --seq1 ${_arg_disassemble_c} \
	    	      --seq2 ${_ptdna} \
		          --out ${_ptdir} \
		          2>$_polap_output_dest"
					done
					shopt -u nullglob # Restore default behavior

					find_folder_with_max_coverage() {
						local max_value=-1
						local max_folder=""

						local folder

						shopt -s nullglob # Enable nullglob so that the loop gets no elements if no match
						for _ptdna in "${_var_mtdna}"/circular_path_*_concatenated.fa; do
							j=$(_polap_utility_extract_number_circular_path "${_ptdna}")
							if [[ -f "${_var_mtdna}/${j}/coverage1.txt" ]]; then
								local value=$(<"${_var_mtdna}/${j}/coverage1.txt")
								if (($(echo "$value > $max_value" | bc -l))); then
									max_value=$value
									max_folder=$j
								fi
							fi
						done
						shopt -u nullglob # Restore default behavior

						if [[ -n $max_folder ]]; then
							echo "$max_folder" # Return the folder number with the largest value
						else
							_polap_log2 "  no valid coverage1.txt files found."
							return 1 # Indicate an error
						fi
					}

					_folder_with_max_coverage=$(find_folder_with_max_coverage)
					if [[ $? -eq 0 ]]; then
						# Calculate absolute difference using bc
						_coverage_ref=$(<"${_var_mtdna}/${_folder_with_max_coverage}/coverage1.txt")
						_coverage_target=$(<"${_var_mtdna}/${_folder_with_max_coverage}/coverage2.txt")
						local diff=$(echo "scale=5; a=(${_coverage_ref} - ${_coverage_target}); if (a < 0) a=-a; a" | bc)

						# Check if difference is less than 0.05
						if (($(echo "$diff < 0.05" | bc -l))); then
							_polap_log2 "  coverage_ref and coverage_target: difference ($diff) is less than 0.05"
							_j_candidate="${_folder_with_max_coverage}"
						fi
					fi

					# select one ptDNA
					if [[ "${_j_candidate}" -ne -1 ]]; then

						cp -p "${_var_mtdna}/${_j_candidate}/coverage1.txt" \
							"${_var_mtdna}/measure_ptdna-blast-alignment-coverage-rate1.txt"
						cp -p "${_var_mtdna}/${_j_candidate}/coverage2.txt" \
							"${_var_mtdna}/measure_ptdna-blast-alignment-coverage-rate2.txt"
						_restarted_fasta="${_var_mtdna}/${_j_candidate}/seq2_restarted.fasta"

						# polish the ptDNA if necessary
						#
						_polap_log3_cmd rm -f "${_arg_unpolished_fasta}"
						_polap_log3_cmd rm -f "${_arg_final_assembly}"

						if [[ -s "${_restarted_fasta}" ]]; then
							_polap_log2 "    ptDNA selected: ${_restarted_fasta}"
							if [[ "${num_circular_paths}" -eq 4 ]] || [[ "${num_circular_paths}" -eq 2 ]]; then
								_polap_log3_cmd cp -p "${_restarted_fasta}" "${_arg_unpolished_fasta}"
							else
								die "    but number of circular patths: ${num_circular_paths} neither 2 nor 4"
							fi
						else
							_polap_log2 "    no such ptDNA: ${_restarted_fasta}"
						fi

						if [[ "${_arg_polish}" == "on" ]]; then
							if [[ -s "${_arg_unpolished_fasta}" ]]; then
								_run_polap_polish
							fi
						fi
					else
						_polap_log2 "    (compare) no ptDNA selected"
					fi
				else
					_polap_log3_cmd cp -p "${_circular_path_fasta}" "${_arg_unpolished_fasta}"
					_polap_log2 "    (inference) ptDNA: ${_circular_path_fasta}"
				fi
			else
				_polap_log2 "  no ptDNA: circular_path count is neither 2 nor 4."
			fi
		fi

		# ptDNA BLAST global alignment
		# pident
		pident=-1
		if _polap_contains_step 15 "${_step_array[@]}"; then
			_polap_log1 "  step 15: percent identity comparison"
			_polap_log2 "    output: ${_var_mtdna}/b/pident.txt"
			_polap_log2 "    output: ${_var_mtdna}/b/blast_results_2.txt"
			_ptdna="${_var_mtdna}/ptdna.1.fa"
			_ptdir="${_var_mtdna}/b"
			if [[ -s "${_var_mtdna}/b/pident.txt" ]] &&
				[[ "${_arg_test}" == "on" ]]; then
				_polap_log2 "  found: ${_var_mtdna}/b/pident.txt, so skipping percent identity ..."
				pident=$(<"${_var_mtdna}/b/pident.txt")
			else
				pident=0
				if [[ -s "${_ptdna}" ]]; then
					_polap_log3_pipe "python $script_dir/run-polap-py-compare2ptdna.py \
    		    --seq1 ${_arg_disassemble_c} \
	    	    --seq2 ${_ptdna} \
		        --out ${_ptdir} \
		        2>$_polap_output_dest"
					if [[ -s "${_var_mtdna}/b/pident.txt" ]]; then
						pident=$(<"${_var_mtdna}/b/pident.txt")
					fi
				else
					_polap_log2 "  no such polished ptDNA assembled: ${_ptdna}"
				fi
			fi
		fi

		# ptDNA summary
		#
		length=-1
		gc=-1
		if _polap_contains_step 16 "${_step_array[@]}"; then
			_polap_log1 "  step 16: ptDNA summory"
			ptdna_file="${_var_mtdna}/ptdna.1.fa"
			if [[ ! -s "${ptdna_file}" ]]; then
				ptdna_file="${_var_mtdna}/ptdna.0.fa"
			fi
			_polap_log2 "    input1: ${ptdna_file}"

			# cflye other steps
			if [[ -s "${ptdna_file}" ]]; then
				# NOTE: subset gfa with the seeds
				# TODO: a circular path

				nseq=$(seqkit stats -Ta "${ptdna_file}" |
					csvtk cut -t -f num_seqs |
					csvtk del-header |
					head -n 1)
				[[ "${nseq}" == 1 ]] || die "ERROR: ${nseq} only one sequence in ${ptdna_file}"
				gc=$(seqkit stats -Ta "${ptdna_file}" |
					csvtk cut -t -f "GC(%)" |
					csvtk del-header |
					head -n 1)
				length=$(seqkit stats -Ta "${ptdna_file}" |
					csvtk cut -t -f sum_len |
					csvtk del-header |
					head -n 1)
			else
				_polap_log2 "  no ptDNA assembled: ${ptdna_file}"
			fi
		fi

		# Append the result to the output file
		# printf "index\tlong\trate\tsize\talpha\tmemory\tnsegments\tnlink\ttotalsegmentlength\tlength\tgc\tcoverage\n" >"${_summary1_file}"
		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
			"${i}" \
			"${sampling_datasize}" \
			"${long_total}" "${long_rate_sample}" "${long_sample_seed}" \
			"${short_total}" "${short_rate_sample}" "${short_sample_seed}" \
			"${_expected_genome_size}" \
			"${peak_ram_size_gb}" \
			"${_pre_alpha}" \
			"${gfa_number_segments}" "${gfa_number_links}" "${gfa_total_segment_length}" \
			"${num_circular_paths}" "${num_circular_nodes}" \
			"$(_polap_get_elapsed_time)" \
			"${_coverage_ref}" \
			"${_coverage_target}" \
			"$pident" \
			"$length" \
			"$gc" \
			"${_j_candidate}" \
			>>"$_summary1_file"

		# Progress
		# Calculate elapsed time and remaining time
		time_per_iteration=$(($(date +%s) - _polap_var_start_time))
		remaining_iterations=$((_total_iterations - i - 1))
		remaining_time=$((remaining_iterations * time_per_iteration))
		status="    iteration $((i + 1))/$_total_iterations, remaining time: $(_polap_get_time_format ${remaining_time})"

		# Display the progress and remaining time on the same line
		if [[ -z "$resume" ]]; then
			if [ "${_arg_verbose}" -eq "1" ]; then
				printf "\r%-${_terminal_width}s" " " >&3
			fi
			_polap_log0_ne "\r$status"
		fi

		# determine the sample size and the Flye's alpha
		_polap_log1 "  end: index $i: sample size: $_sampling_datasize_bp alpha: $_pre_alpha memory: ${peak_ram}"
	done

	if [[ -z "$resume" ]]; then
		_polap_log0 ""
	fi

	if [[ "${_arg_test}" == "off" ]]; then
		rm -f "${_short_read1}"
		rm -f "${_long_read}"
	fi

	return 0
}
