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

# arg1: input fasta or fastq file
# arg2: output text file with the total number of bases in the input sequence
# file
_polap_lib_fastq-total-length-of() {
	local input1="$1"
	local output="$2"

	check_file_existence "${input1}"

	_polap_log3_pipe "seqkit stats -Ta ${input1} -j ${_arg_threads} |\
		csvtk cut -t -f sum_len |\
		csvtk del-header \
		>${output}"

	return 0
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
_polap_lib_fastq-concatenate-fastq-files() {
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
			_polap_log2 "    concatenated fastq file: $input_file1 + $input_file2 -> $output_file"
		else
			_polap_log2 "    concatenated fastq file: $input_file1 -> $output_file"
		fi
	)

}
