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

# Function to determine the minimum of two values
# Each of the two files has a single line of numbers.
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
	local i

	if ((count == 1)); then
		# If count is 1, create an array with a single item using the start value
		array=("$start")
	else
		# Calculate the step size
		local step=$(((end - start) / (count - 1)))

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
# rstatus="$?"
# if [[ "$rstatus" -eq 0 ]]; then
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

# Estimate the genome size from a short-read dataset
# input1: a short-read fastq data
# input2: an output folder
# output: the genome size
_polap_find-genome-size() {
	local _short_read1=$1
	local _outdir=$2
	local rstatus

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
	local _outdir_genome_size="${_outdir}/short_expected_genome_size.txt"
	local _outdir_jellyfish_out="${_outdir}/jellyfish_out"
	local _outdir_jellyfish_out_histo="${_outdir}/jellyfish_out.histo"

	# See https://bioinformatics.uconn.edu/genome-size-estimation-tutorial/
	if [ -s "${_outdir_genome_size}" ] && [ "${_arg_redo}" = "off" ]; then
		_polap_log2 "  found: ${_outdir_genome_size}, so skipping the genome size estimation ..."
		_polap_log3_file "${_outdir_genome_size}"
	else
		if [ -s "${_outdir_jellyfish_out}" ] && [ "${_arg_redo}" = "off" ]; then
			_polap_log2 "  found: ${_outdir_jellyfish_out}, so skipping the JellyFish counting ..."
			_polap_log3_file "${_outdir_jellyfish_out}"
		else
			local unzipped_file
			unzipped_file=$(_polap_gunzip_file "${_short_read1}")
			rstatus="$?"
			if [[ "$rstatus" -eq 0 ]]; then
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

################################################################################
# It is not as straightforward as it might seem. The simplest approach is to test plastid genome assembly by subsampling long-read data. A more complex approach involves using a slightly modified version of Flye, called cflye, which provides coverage options for generating disjointigs. Short-read data is still used to get a rough estimate of the overall genome size, but once subsampling is applied, the size of the genome assembled from long reads can change. The sizes of both the long-read and short-read data are defined by their total number of bases. When subsampling, we reduce the amount of short-read data proportionally to match the size of the subsampled long-read data, then use JellyFish (a k-mer counter) to predict the total size of the genome assembled from these long reads.
#
# If the total long-read data is smaller than the short-read data, we simply subsample the short reads to match the long-read dataset size. However, if the long-read data is larger than the short-read data, we use all the short reads. Because our goal is to estimate the approximate genome size that can be assembled from subsampled long reads, some level of error is to be expected.
#
# In the Flye genome assembly program, disjointigs are generated, and those whose read coverage falls below a certain cutoff are filtered out. We set this initial cutoff using the size of the long-read data (A = 30 Mb, where Mb stands for megabases) and alpha (α = 0.1). We then compare the total length (S) of the disjointigs against a fixed value M. If S is greater than M, alpha is increased by an amount delta (δ = 0.1); if S is less than M, alpha is decreased by the same delta. Along with this process, after each cycle we also increase the long-read data size A by d (for example, 10 Mb) and repeat these cycles until A reaches B, for a total of n cycles. The values of A, B, and n are predetermined, and d is derived from these values. For instance, if A = 30 Mb, B = 120 Mb, and n = 10, then d = 10 Mb.
#
# At the end of each cycle, we extract and evaluate potential plastid genome candidates from the newly assembled genome. To qualify as a plastid candidate, the assembly must consist of a single circular contig or three contigs that form a circular structure. We then assess the size of this circular contig, the number of expected genes, and so on. These values can be plotted as Y-coordinates against the cycle number on the X-axis in a scatter plot, allowing us to visualize progress across cycles.
################################################################################
function _run_polap_disassemble {
	# Enable debugging if DEBUG is set
	[ "$DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_/Menu_/)"

	source "$script_dir/polap-variables-common.sh"
	local _disassemble_dir="${_arg_outdir}/disassemble"
	local _summary_table="${_disassemble_dir}/summary.txt"

	help_message=$(
		cat <<HEREDOC
Plastid genome assembly by subsampling long-read data without reference

Inputs
------

- long-read data: ${_arg_long_reads} (default: l.fq)
- short-read data 1 (optional): ${_arg_short_read1} (no default)
- short-read data 2 (optional): ${_arg_short_read2} (no default)

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
--disassemble-n ${_arg_disassemble_n}: the number of steps (max 1000)
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

--genomesize <NUMBER>: the number of steps is equal to 1.
--genomesize-a ${_arg_genomesize_a}: the smallest genome-size step size
--genomesize-b ${_arg_genomesize_b}: the largest genome-size step size
--genomesize-n ${_arg_genomesize_n}: the number of steps in the genome-size

Stage 2:
--disassemble-s: the threshold to change alpha
--disassemble-alpha ${_arg_disassemble_alpha}: the minimum Flye's disjointig coverage

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
cp o/00-bioproject/2-mtdna.fasta ptdna-epauciflora.fa
2:
Stage 1 - long-read and short-read data
$(basename "$0") disassemble -l SRR7153095.fastq -a SRR7161123_1.fastq -b SRR7161123_2.fastq --disassemble-compare-to-fasta ptdna-epauciflora.fa --disassemble-min-memory 9
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

HEREDOC
	)

	# Display help message
	[[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]] && _polap_echo0 "${help_message}" && return 0
	[[ ${_arg_menu[1]} == "redo" ]] && _arg_redo="on"

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		case "${_arg_menu[2]}" in
		'' | *[!0-9]*)
			_polap_log0_column "${_summary_table}"
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

	local rstatus
	local i # important in Bash scripting
	_arg_plastid="on"
	_arg_disassemble_s=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_s}")
	_arg_disassemble_a=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_a}")
	_arg_disassemble_b=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_b}")
	_arg_disassemble_m=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_m}")
	_arg_disassemble_s_max=$(_polap_utility_convert_unit_to_bp "${_arg_disassemble_s_max}")

	if [[ -n "${_arg_genomesize}" ]]; then
		_arg_genomesize=$(_polap_utility_convert_unit_to_bp "${_arg_genomesize}")
	fi
	local _input_short_reads="${_arg_outdir}/s1.fq"

	_polap_log0 "Assemble plastid genomes by subsampling long-read data"
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

	# NOTE: the processing of short-read data used to be in step-disassemble.
	# because the repetition of step-disassemble, we process the short-read data
	# here in this function caller.
	#
	# make a short-read dataset if there are short-read data
	# Add commands for step 2
	if [[ "${_arg_short_read1_is}" == "on" ]] ||
		[[ "${_arg_short_read2_is}" == "on" ]]; then
		(
			_polap_log1 "  prestep 1: concatenate the short-read data"
			# create o/s1.fq to deal with a single or two short-read files
			# overwrite the s1.fq
			prestep1() {
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
			prestep1
		)

		_polap_log1 "  prestep 2: count the total number of bases in the short-read data"
		# Add commands for step 3
		(
			prestep2() {
				# step 3
				if [ -s "${_polap_var_outdir_short_total_length}" ] && [ "${_arg_redo}" = "off" ]; then
					_polap_log2 "  found: ${_polap_var_outdir_short_total_length}, skipping total-length-short ..."
				else
					_polap_total-length-short \
						"${_input_short_reads}" \
						"${_polap_var_outdir_short_total_length}"

					local _total_length=$(<"${_polap_var_outdir_short_total_length}")
					local _total_short_read=$(_polap_utility_convert_bp ${_total_length})
					_polap_log1 "  total length of the short-read dataset (bases): ${_total_short_read}"
				fi
			}
			prestep2
		)
	fi

	# after preparation of short-read datasets
	if [[ -z "${_arg_disassemble_s}" ]]; then
		_polap_log1 "  stage 1: assembling with a range of many parameters"
		_polap_log1 "  no option: --disassemble-s"
		_polap_log1 "    -> use a range of rate for subsampling long-read sequencing data"
		if [[ "${_arg_short_read1_is}" == "off" ]] &&
			[[ "${_arg_short_read2_is}" == "off" ]]; then
			_polap_log1 "  no short-read data, so use a range of genome sizes"
			if [[ -z "${_arg_steps_include}" ]]; then
				_arg_steps_include="1-15"
				_arg_steps_exclude="2,3,6"
			fi

			# NOTE: ptDNA varation over the fixed genome size
			# genome sizes: if we cannot estimate the genome size, then use a list
			# of variable genome sizes in the alpha iteration
			# for each fix genome size, we check the distribution of genome lengths
			# over the sampling size. Then, we choose one genome size that shows
			# a better behavior in the ...
			# we could check the variance of the ptDNA's features depending on the
			# fixed genome size.
			if [[ -n "${_arg_genomesize}" ]]; then
				_arg_genomesize_a="${_arg_genomesize}"
				_arg_genomesize_b="${_arg_genomesize}"
				_arg_genomesize_n=1
			fi
			local genomesize_array=($(_polap_array_numbers_between_two \
				"$_arg_genomesize_a" \
				"$_arg_genomesize_b" \
				"$_arg_genomesize_n"))
			local genomesize_array_length=${#genomesize_array[@]}
			_polap_log1 "  genome sizes (${_arg_genomesize_n}): ${_arg_genomesize_a} ~ ${_arg_genomesize_b}"

			_polap_log3_cmd mkdir -p "${_disassemble_dir}"
			_polap_log3_cmd rm -f "${_disassemble_dir}/index-genomesize.txt"
			# local i=0 this does not work. It must be local
			# at the beginning of the function.
			for ((i = 0; i < genomesize_array_length; i++)); do
				local genomesize="${genomesize_array[$i]}"
				printf "%d\t%d\n" $i $genomesize \
					>>"${_arg_outdir}/disassemble/index-genomesize.txt"
			done

			# local i=0 this does not work. It must be local
			# at the beginning of the function.
			for ((i = 0; i < genomesize_array_length; i++)); do
				local genomesize="${genomesize_array[$i]}"
				_run_polap_step-disassemble "${i}" "${genomesize}"
				rstatus="$?"
				_polap_log0 "genomesize: rstatus: $rstatus"
				_polap_log0 "genomesize: i: $i"
				echo "genomesize: rstatus: $rstatus" >&2
				echo "genomesize: i: $i" >&2
				[[ "$rstatus" -ne 0 ]] && return "$rstatus"
			done

			# concatenate all summary1.txt to the summary table.
			i=0
			awk 'NR == 1' "${_disassemble_dir}/${i}/summary1.txt" \
				>"${_summary_table}"
			for ((i = 0; i < genomesize_array_length; i++)); do
				awk 'NR > 1' "${_disassemble_dir}/${i}/summary1.txt" \
					>>"${_summary_table}"
			done
		else
			_polap_log1 "  genome size estimated by the short-read data"
			if [[ -z "${_arg_steps_include}" ]]; then
				_arg_steps_include="1-15"
				_arg_steps_exclude="2,3"
			fi
			check_file_existence "${_input_short_reads}"
			_run_polap_step-disassemble 0
			rstatus="$?"
			[[ "$rstatus" -ne 0 ]] && return "$rstatus"
		fi
	else
		_polap_log1 "  stage 2: assembling with a single parameter"
		local assembled_fasta="${_disassemble_dir}/x/0/ptdna.fa"
		_arg_unpolished_fasta="pt.0.fasta"
		_arg_final_assembly="pt.1.fa"
		_arg_contigger="off"
		_arg_disassemble_a="${_arg_disassemble_s}"
		_arg_disassemble_b="${_arg_disassemble_s}"
		_arg_disassemble_b_is="on"
		_arg_disassemble_n=1
		_polap_log1 "  input1: --disassemble-s=${_arg_disassemble_s}"
		if [[ "${_arg_short_read1_is}" == "off" ]] &&
			[[ "${_arg_short_read2_is}" == "off" ]]; then
			if [[ -z "${_arg_steps_include}" ]]; then
				_arg_steps_include="1-15"
				_arg_steps_exclude="2,3,6"
			fi
			if [[ -z "${_arg_genomesize}" ]]; then
				_polap_log0 "ERROR: menu disassemble without short-read data and with --disassemble-s option"
				return "${_POLAP_ERR_CMD_OPTION_GENOMESIZE}"
			else
				_polap_log1 "  input2: genome size: ${_arg_genomesize}"
			fi

			_run_polap_step-disassemble
			rstatus="$?"
			[[ "$rstatus" -ne 0 ]] && return "$rstatus"

			# percentile in the length distribution
			# need --anotherdir
			_run_polap_disassemble-stats
			_polap_log3_cmd cp "${assembled_fasta}" "${_arg_final_assembly}"
		else
			if [[ -z "${_arg_steps_include}" ]]; then
				_arg_steps_include="1-15"
				_arg_steps_exclude="2,3"
			fi
			_polap_log1 "  genome size estimated by the short-read data"
			check_file_existence "${_input_short_reads}"
			_run_polap_step-disassemble
			rstatus="$?"
			[[ "$rstatus" -ne 0 ]] && return "$rstatus"

			# percentile in the length distribution
			# need --anotherdir
			_run_polap_disassemble-stats
			_run_polap_prepare-polishing

			_polap_log3_cmd cp "${assembled_fasta}" "${_arg_unpolished_fasta}"
			_run_polap_polish
		fi
		_polap_log0 "Final plastid genome sequence: ${_arg_final_assembly}"
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

	_arg_plastid="on"
	local i=0
	local genomesize

	local _input_short_reads="${_arg_outdir}/s1.fq"

	if [[ $# -gt 0 ]]; then
		local _disassemble_dir="${_arg_outdir}/disassemble/${1}"
	else
		local _disassemble_dir="${_arg_outdir}/disassemble/x"
	fi

	if [[ $# -gt 1 ]]; then
		genomesize="${2}"
	else
		genomesize=""
		if [[ -n "${_arg_genomesize}" ]]; then
			genomesize="${_arg_genomesize}"
		fi
	fi

	_polap_log0 "  assemble a plastid genome without reference"
	_polap_log1 "  input1: long-read: ${_arg_long_reads}"
	if [[ "${_arg_short_read1_is}" == "off" ]]; then
		_polap_log1 "  input2: no short-read1"
	else
		_polap_log1 "  input2: short-read1: ${_arg_short_read1}"
	fi
	if [[ "${_arg_short_read2_is}" == "off" ]]; then
		_polap_log1 "  input3: no short-read2"
	else
		_polap_log1 "  input3: short-read2: ${_arg_short_read2}"
	fi
	_polap_log1 "  input4: option --plastid: ${_arg_plastid}"
	if [[ -n "${genomesize}" ]]; then
		_polap_log1 "  input5: genomesize: ${genomesize}"
	else
		_polap_log1 "  input5: no given genomesize"
	fi

	if [[ -n "${_arg_species}" ]]; then
		_polap_log1 "  input6: species: ${_arg_species}"
		_run_polap_get-mtdna
		_arg_disassemble_compare_to_fasta="${_arg_outdir}/00-bioproject/2-mtdna.fasta"
	elif [[ -n "${_arg_disassemble_compare_to_fasta}" ]]; then
		_polap_log1 "  input6: compare to fasta: ${_arg_disassemble_compare_to_fasta}"
		check_file_existence "${_arg_disassemble_compare_to_fasta}"
	elif [[ -s "${_polap_var_project_mtdna_fasta2}" ]]; then
		local n=$(grep ">" "${_polap_var_project_mtdna_fasta2}" | wc -l)
		if [[ "${n}" == "1" ]]; then
			_arg_disassemble_compare_to_fasta="${_polap_var_project_mtdna_fasta2}"
			_polap_log1 "  input6: compare to fasta: ${_arg_disassemble_compare_to_fasta}"
		else
			_polap_log1 "  input6: neither species nor fasta file is given"
		fi
	else
		_polap_log1 "  input6: neither species nor fasta file is given"
		_polap_log1 "    one of the sequences of an assembly graph is chosen at random"
	fi
	_polap_log1 "  output: ${_disassemble_dir}"

	# select steps
	local include="${_arg_steps_include}"
	local exclude="${_arg_steps_exclude}" # Optional range or list of steps to exclude
	local step_array=()
	local exclude_array=()

	step_array=($(_polap_parse_steps "${include}" "${exclude}"))
	if [[ "$?" -ne 0 ]]; then
		_polap_log0 "ERROR: Error parsing steps."
		return "${_POLAP_ERR_CMD_OPTION_STEPS}"
	fi

	_polap_log1 "  executing steps: ${step_array[*]}"

	if _polap_contains_step 1 "${step_array[@]}"; then
		_polap_log1 "  step 1: count the total number of bases in the long-read data"
		(
			step1() {
				# step 1
				if [ -s "${_polap_var_outdir_long_total_length}" ] && [ "${_arg_redo}" = "off" ]; then
					_polap_log2 "  found: ${_polap_var_outdir_long_total_length}, skipping total-length-long ..."
				else
					_run_polap_total-length-long
				fi
			}
			step1
		)
	fi

	# NOTE: it must have been done in the function caller
	if _polap_contains_step 2 "${step_array[@]}"; then
		_polap_log1 "  no step 2: concatenate the short-read data"
		check_file_existence "${_input_short_reads}"
	fi

	if _polap_contains_step 3 "${step_array[@]}"; then
		_polap_log1 "  no step 3: count the total number of bases in the short-read data"
		check_file_existence "${_input_short_reads}"
	fi

	if _polap_contains_step 4 "${step_array[@]}"; then
		_polap_log1 "  step 4: prepare for the loop iteration"
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
		_polap_log1 "    the disassemble A (bp): ${_arg_disassemble_a}"
		_polap_log1 "    the disassemble B (bp): ${_arg_disassemble_b}"
		_polap_log1 "    the disassemble N: ${_arg_disassemble_n}"
		_polap_log1 "    the disassemble M (bp): ${_arg_disassemble_m}"
		_polap_log1 "    the disassemble alpha: ${_arg_disassemble_alpha}"
		_polap_log1 "    the disassemble delta: ${_arg_disassemble_delta}"
		if ((_arg_disassemble_a <= _arg_disassemble_b)); then
			_polap_log2 "  check: ${_arg_disassemble_a} is less than or equal to ${_arg_disassemble_b}."
		else
			die "ERROR: ${_arg_disassemble_a} (A) is greater than ${_arg_disassemble_b} (B)."
		fi

		datasize_array=($(_polap_array_numbers_between_two \
			"$_arg_disassemble_a" \
			"$_arg_disassemble_b" \
			"$_arg_disassemble_n"))
		_polap_log2 "  generated step array (${#datasize_array[@]}): ${datasize_array[@]:0:3} ... ${datasize_array[@]: -2}"

		if [ -d "${_disassemble_dir}" ] && [ "${_arg_redo}" = "on" ]; then
			_polap_log3_cmd rm -rf "${_disassemble_dir}"
		fi
		_polap_log3_cmd mkdir -p "${_disassemble_dir}"

		datasize_array_length=${#datasize_array[@]}

		rm -f "${_disassemble_dir}/index.txt"
		for ((i = 0; i < datasize_array_length; i++)); do
			sampling_datasize="${datasize_array[$i]}"
			printf "%d\t%d\n" $i $sampling_datasize >>"${_disassemble_dir}/index.txt"
		done

		local _alpha="${_arg_disassemble_alpha}"

		# Iterate using indices of the data size
		local _size_index="${datasize_array_length}"
		if ((_arg_end_index > 0)); then
			_size_index="${_arg_end_index}"
		fi
	fi

	local _short_read1=${_arg_outdir}/s.fq
	local _long_read=${_arg_outdir}/l.fq.gz
	local is_stop="off"
	local summary1_file="${_disassemble_dir}/summary1.txt"
	local num_circular_paths
	local num_circular_nodes
	if _polap_contains_step 5 "${step_array[@]}"; then
		# Add header to the output file
		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
			"index" \
			"long" \
			"rate" \
			"size" \
			"alpha" \
			"genomesize" \
			"memory" \
			"nsegments" \
			"nlink" \
			"totalsegment" \
			"ncircle" \
			"ncirclenode" \
			"length" \
			"gc" \
			"coverage" \
			>"${summary1_file}"
		# printf "index\trate\tsize\talpha\tmemory\tnsegments\tnlink\ttotalsegmentlength\n" >"${summary1_file}"
	else
		is_stop="on"
	fi
	local i_last=0
	for ((i = _arg_start_index; i < _size_index; i++)); do
		if [[ "${is_stop}" == "on" ]]; then
			break
		fi
		i_last=$((i_last + 1))
		local sampling_datasize="${datasize_array[$i]}"
		local preAlpha="${_alpha}"
		local _outdir=${_disassemble_dir}/$i
		local _wga=${_outdir}

		local _sampling_datasize_bp=$(_polap_utility_convert_bp ${sampling_datasize})
		_polap_log1 "($i) sample data size: ${_sampling_datasize_bp}"
		_polap_log1 "  alpha: ${_alpha}"
		_polap_log3_cmd mkdir -p "${_outdir}"

		if _polap_contains_step 5 "${step_array[@]}"; then
			_polap_log1 "  step 5: subsample long-read data"
			# subsample the long-read data
			local long_total=$(cat "${_polap_var_outdir_long_total_length}")

			_polap_log2 "  output: ${_long_read}"
			RATE=$(echo "scale=10; ${sampling_datasize}/${long_total}" | bc)
			local rate_sample_long="${RATE}"
			local seed=${_arg_random_seed:-$RANDOM}
			_polap_log2 "    random seed: ${seed}"
			_polap_log3_pipe "seqkit sample \
		      -p ${RATE} \
		      -s ${seed} \
		      ${_arg_long_reads} \
		      -o ${_long_read} 2>${_polap_output_dest}"
		fi

		if _polap_contains_step 6 "${step_array[@]}"; then
			_polap_log1 "  step 6: subsample the short-read data: ${sampling_datasize} (bp)"
			check_file_existence "${_input_short_reads}"

			local short_total=$(cat "${_polap_var_outdir_short_total_length}")

			RATE=$(echo "scale=10; ${sampling_datasize}/${short_total}" | bc)
			local seed=${_arg_random_seed:-$RANDOM}
			_polap_log2 "    random seed: ${seed}"
			_polap_log3_pipe "seqkit sample \
        -p ${RATE} \
        -s ${seed} \
        ${_input_short_reads} \
        -o ${_short_read1} 2>${_polap_output_dest}"
		fi

		if _polap_contains_step 7 "${step_array[@]}"; then
			_polap_log1 "  step 7: use a given genome size or estimate one using the subsample of the short-read data"
			local _expected_genome_size

			if [[ "${_arg_short_read1_is}" == "off" ]] &&
				[[ "${_arg_short_read2_is}" == "off" ]]; then
				_expected_genome_size="${genomesize}"
				_polap_log2 "  given genome size: ${_expected_genome_size} bp"
			else
				_polap_find-genome-size \
					"${_short_read1}" \
					"${_outdir}"
				_expected_genome_size=$(cat "${_outdir}/short_expected_genome_size.txt")
				_polap_log2 "  estimated genome size: ${_expected_genome_size} bp"
			fi
		fi

		if _polap_contains_step 8 "${step_array[@]}"; then
			_polap_log1 "  step 8: cFlye on the subsample of the long-read data"
			# cflye assemble stage
			_polap_log1 "  executing the whole-genome assembly using cflye ... be patient!"
			local _command1="/home/goshng/all/polap/Flye/bin/cflye \
          ${_arg_flye_data_type} \
          ${_long_read} \
			    --out-dir ${_wga} \
			    --disjointig-min-coverage ${_alpha} \
			    --threads ${_arg_threads}"

			_command1+=" \
			  --asm-coverage ${_arg_flye_asm_coverage} \
			  --genome-size ${_expected_genome_size}"
			if [[ "${_arg_flye_asm_coverage}" -le 0 ]]; then
				die "ERROR: must be positive --flye-asm-coverage ${_arg_flye_asm_coverage}"
			fi

			_command1+=" \
		    --stop-after contigger"

			if [[ "${_arg_debug}" = "on" ]]; then
				_command1+=" \
		      --debug"
			fi

			_command1+=" \
		    2>${_polap_output_dest}"

			if [[ -d "${_outdir}/00-assembly" ]]; then
				_polap_log2 "  found: ${_outdir}/00-assembly/draft_assembly.fasta, ..."
			else
				_polap_log3_pipe "${_command1}"
			fi

			# use all the data in the long-read not just the sampled one
			# for the polishing stage
			if [[ "${_arg_contigger}" == "off" ]] &&
				[[ -d "${_outdir}/30-contigger" ]]; then
				local _command1="/home/goshng/all/polap/Flye/bin/cflye \
          ${_arg_flye_data_type} \
          ${_arg_long_reads} \
			    --out-dir ${_wga} \
			    --resume \
			    --threads ${_arg_threads}"

				if [[ "${_arg_debug}" = "on" ]]; then
					_command1+=" \
		        --debug"
				fi

				_command1+=" \
		      2>${_polap_output_dest}"

				_polap_log3_pipe "${_command1}"
			else
				_polap_log2 "  conttiger is on or not found: ${_outdir}/30-contigger"
			fi

			local peak_ram_size_gb
			local peak_ram
			if [[ -d "${_outdir}/30-contigger" ]]; then
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
				if [[ -n "${_arg_disassemble_s_max}" ]]; then
					if [[ "${sampling_datasize}" -gt "${_arg_disassemble_s_max}" ]]; then
						is_stop="on"
					fi
				fi
			else
				_polap_log2 "  cFlye: no genome assembly"
				peak_ram_size_gb=-1
				peak_ram=-1
			fi
		fi

		local _contigger_edges_gfa="${_outdir}/30-contigger/graph_final.gfa"
		if [[ "${_arg_contigger}" == "off" ]]; then
			if [[ -s "${_outdir}/assembly_graph.gfa" ]]; then
				_polap_log2 "  backup ${_contigger_edges_gfa}"
				_polap_log3_cmd cp -p "${_contigger_edges_gfa}" "${_contigger_edges_gfa}.backup"
				_polap_log2 "  use ${_outdir}/assembly_graph.gfa by copying the contigger graph_final.gfa"
				_polap_log3_cmd cp -p "${_outdir}/assembly_graph.gfa" "${_contigger_edges_gfa}"
			fi
		fi

		if _polap_contains_step 9 "${step_array[@]}"; then
			_polap_log1 "  step 9: adjust the alpha"
			# adjust alpha
			if [[ -s "${_outdir}/00-assembly/draft_assembly.fasta" ]]; then
				local _outdir_draft_assembly_size="${_outdir}/draft_assembly_size.txt"
				_polap_log3_pipe "seqkit stats -Ta ${_outdir}/00-assembly/draft_assembly.fasta |\
		      csvtk cut -t -f sum_len |\
		      csvtk del-header \
		      >${_outdir_draft_assembly_size}"

				local _draft_assembly_size=$(cat "${_outdir_draft_assembly_size}")
				local _draft_assembly_size_bp=$(_polap_utility_convert_bp ${_draft_assembly_size})
				_polap_log1 "  assembly for data size: ${_draft_assembly_size_bp}"
				local _arg_disassemble_m_bp=$(_polap_utility_convert_bp ${_arg_disassemble_m})

				if ((_arg_disassemble_m < _draft_assembly_size)); then
					_alpha=$(_polap_adjust_alpha_by_delta "${_alpha}" \
						"${_arg_disassemble_delta}" \
						"increase")
					_polap_log1 "  Alpha increased: ${_alpha}"
				else
					_alpha=$(_polap_adjust_alpha_by_delta "${_alpha}" \
						"${_arg_disassemble_delta}" \
						"decrease")
					_polap_log1 "  Alpha decrease (if not negative): ${_alpha}"
				fi
			else
				local _sampling_datasize_bp=$(_polap_utility_convert_bp ${sampling_datasize})
				_polap_log1 "  no assembly for data size: ${_sampling_datasize_bp}"

				_alpha=$(_polap_adjust_alpha_by_delta "${_alpha}" \
					"${_arg_disassemble_delta}" \
					"decrease")
				_polap_log1 "  Alpha decrease (if not negative): ${_alpha}"
			fi
		fi

		# cflye other steps
		if [[ -s "${_contigger_edges_gfa}" ]]; then
			if _polap_contains_step 10 "${step_array[@]}"; then
				_polap_log1 "  step 10: summary statistics of gfa: ${_contigger_edges_gfa}"
				gfatools stat -l "${_contigger_edges_gfa}" >"${_outdir}/graph_final.txt" 2>"$_polap_output_dest"
				local num_segments=$(grep "Number of segments:" "${_outdir}/graph_final.txt" | awk '{print $4}')
				local num_links=$(grep "Number of links:" "${_outdir}/graph_final.txt" | awk '{print $4}')
				local total_segment_length=$(grep "Total segment length:" "${_outdir}/graph_final.txt" | awk '{print $4}')
			fi

			if _polap_contains_step 11 "${step_array[@]}"; then
				_polap_log1 "  step 11: annotate"
				# find the contig with more plastid genes than mitochondrial
				local _ga="${_outdir}"
				local _ga_annotation_all="${_ga}/assembly_info_organelle_annotation_count-all.txt"
				polap_annotate "${_contigger_edges_gfa}" "${_ga_annotation_all}"
			fi

			local _ga="${_outdir}"
			local _mtcontigname="${_ga}/mt.contig.name"
			# create mt.contig.name.txt
			if _polap_contains_step 12 "${step_array[@]}"; then
				_polap_log1 "  step 12: seeds"
				_polap_log1 "    output: ${_mtcontigname}"
				# depth-filtering
				polap_disassemble-seeds "${_contigger_edges_gfa}" \
					"${_ga_annotation_all}" \
					"${_mtcontigname}"
			fi

			# NOTE: subset gfa with the seeds
			# TODO: a circular path

			if _polap_contains_step 13 "${step_array[@]}"; then
				_polap_log1 "  step 13: extract ptDNA sequences from the gfa"
				# locate the circular paths
				# extract candidate circular ptDNA sequences
				_polap_log2 "    input1: ${_contigger_edges_gfa}"
				_polap_log2 "    input2: ${_mtcontigname}"
				_polap_log2 "    output: ${_outdir}"
				_polap_log3_pipe "python \
          $script_dir/run-polap-py-find-plastid-gfa2fasta.py \
		      --gfa ${_contigger_edges_gfa} \
		      --seed ${_mtcontigname} \
		      --out ${_outdir} \
		      2>$_polap_output_dest"
				num_circular_paths=$(<"${_outdir}/circular_path_count.txt")
				num_circular_nodes=$(<"${_outdir}/circular_path_nodes.txt")
				_polap_log3_cat "${_outdir}/circular_path.txt"
			else
				num_circular_paths=0
				num_circular_nodes=0
			fi

			local _ptdna="${_outdir}/circular_path_1_concatenated.fa"
			if [[ -n "${_arg_disassemble_compare_to_fasta}" ]] &&
				[[ -s "${_ptdna}" ]] &&
				_polap_contains_step 14 "${step_array[@]}"; then
				_polap_log1 "  step 14: select the best ptDNA using FASTA: ${_arg_disassemble_compare_to_fasta}"

				# for each reference ptDNA,
				#   select one candidate ptDNA
				#   rearrange it so that we could do a pairwise sequence alignment
				_polap_log2 "    input1: ${_arg_disassemble_compare_to_fasta}"
				_polap_log2 "    output: ${_outdir}"
				local j=0
				for j in 1 2 3 4; do
					local _ptdna="${_outdir}/circular_path_${j}_concatenated.fa"
					local _ptdir="${_outdir}/${j}"

					if [[ -s "${_ptdna}" ]]; then
						_polap_log3_pipe "python $script_dir/run-polap-py-compare2ptdna.py \
    		      --seq1 ${_arg_disassemble_compare_to_fasta} \
	    	      --seq2 ${_ptdna} \
		         --out ${_ptdir} \
		         2>$_polap_output_dest"
					else
						_polap_log2 "  no such ptDNA assembled: ${_ptdna}"
					fi
				done

				find_max_folder_with_coverage() {
					local _outdir="$1"
					local max_value=-1
					local max_folder=""

					local folder
					for folder in 1 2 3 4; do
						if [[ -f "${_outdir}/${folder}/coverage.txt" ]]; then
							local value=$(cat "${_outdir}/${folder}/coverage.txt")
							if (($(echo "$value > $max_value" | bc -l))); then
								max_value=$value
								max_folder=$folder
							fi
						# else
						# 	_polap_log2 "  warning: no such file: ${_outdir}/${folder}/coverage.txt"
						fi
					done

					if [[ -n $max_folder ]]; then
						echo "$max_folder" # Return the folder number with the largest value
					else
						_polap_log2 "  no valid coverage.txt files found."
						return 1 # Indicate an error
					fi
				}

				folder_with_max_coverage=$(find_max_folder_with_coverage "${_outdir}")
				if [[ $? -eq 0 ]]; then
					_polap_log2 "  Folder with max coverage: $folder_with_max_coverage"
					_polap_log2 "  sequence: ${_outdir}/${folder_with_max_coverage}/seq2_restarted.fasta"
					_polap_log2_cat "${_outdir}/${folder_with_max_coverage}/coverage.txt"
					cp "${_outdir}/${folder_with_max_coverage}/seq2_restarted.fasta" "${_outdir}/ptdna.fa"
					cp "${_outdir}/${folder_with_max_coverage}/coverage.txt" "${_outdir}/measure_ptdna-blast-alignment-coverage-rate.txt"
				else
					_polap_log2 "  no max coverage."
					local _ptdna="${_outdir}/circular_path_1_concatenated.fa"
					_polap_log1 "  ptDNA selected: ${_ptdna}"
					_polap_log3_cmd cp "${_ptdna}" "${_outdir}/ptdna.fa"
					echo 0 >"${_outdir}/measure_ptdna-blast-alignment-coverage-rate.txt"
				fi
			else
				local _ptdna="${_outdir}/circular_path_1_concatenated.fa"
				_polap_log1 "  step 14: no species is given by --species option"
				_polap_log1 "    no species is given by --species option"
				if [[ -s "${_ptdna}" ]] &&
					[[ "${num_circular_paths}" -eq 4 ]] &&
					[[ "${num_circular_paths}" -eq 2 ]]; then
					_polap_log1 "    ptDNA selected: ${_ptdna}"
					_polap_log3_cmd cp "${_ptdna}" "${_outdir}/ptdna.fa"
				else
					_polap_log1 "  no such output ptDNA: ${_ptdna}"
					_polap_log1 "  no such output ptDNA: ${_outdir}/ptdna.fa"
				fi
			fi

			local ptdna_file="${_outdir}/ptdna.fa"
			if [[ -s "${ptdna_file}" ]]; then
				local nseq=$(seqkit stats -Ta "${ptdna_file}" |
					csvtk cut -t -f num_seqs |
					csvtk del-header |
					head -n 1)
				[[ "${nseq}" == 1 ]] || die "ERROR: ${nseq} only one sequence in ${ptdna_file}"
				local gc=$(seqkit stats -Ta "${ptdna_file}" |
					csvtk cut -t -f "GC(%)" |
					csvtk del-header |
					head -n 1)
				local length=$(seqkit stats -Ta "${ptdna_file}" |
					csvtk cut -t -f sum_len |
					csvtk del-header |
					head -n 1)
			else
				local length=0
				local gc=0
			fi

			if [[ -s "${_outdir}/measure_ptdna-blast-alignment-coverage-rate.txt" ]]; then
				local coverage=$(cat "${_outdir}/measure_ptdna-blast-alignment-coverage-rate.txt")
			else
				local coverage=0
			fi

		else
			local num_segments=0
			local num_links=0
			local total_segment_length=0
			local length=0
			local gc=0
			local coverage=0
		fi

		# Append the result to the output file
		# printf "index\tlong\trate\tsize\talpha\tmemory\tnsegments\tnlink\ttotalsegmentlength\tlength\tgc\tcoverage\n" >"${summary1_file}"
		printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
			"$i" "${long_total}" "${rate_sample_long}" "$sampling_datasize" \
			"$preAlpha" "${_expected_genome_size}" "${peak_ram_size_gb}" \
			"$num_segments" "$num_links" "$total_segment_length" \
			"$num_circular_paths" "$num_circular_nodes" \
			"$length" "$gc" "$coverage" \
			>>"$summary1_file"

		# determine the sample size and the Flye's alpha
		_polap_log1 "status: index $i: length: $_sampling_datasize_bp alpha: $preAlpha ram: ${peak_ram}"
	done

	return 0
}
