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

# Function to convert a comma-delimited string to a Bash array
# input="A_1.fastq,A_2.fastq,A_3.fastq"
# result=($(convert_elements_to_array "$input"))
_polap_lib_convert-list-to-array() {
	local input_string=$1
	IFS=',' read -r -a array <<<"$input_string"
	echo "${array[@]}"
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
_polap_lib_array-numbers-a2b() {
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

# function _disassemble-make-index-for-p
#
# arg1: index.table (output)
# arg2: number of the sampling steps
# arg4: min sample size
# arg3: max sample size
#
# index for short-read p
function _polap_lib_array-make-index {
	local _index_table="${1}"
	local _disassemble_n="${2}"
	local _disassemble_a="${3}"
	local _disassemble_b="${4}"

	_polap_log3 "    create a loop index: ${_index_table}"
	_polap_log3 "      the disassemble N: ${_disassemble_n}"
	_polap_log3 "      the disassemble A (bp): ${_disassemble_a}"
	_polap_log3 "      the disassemble B (bp): ${_disassemble_b}"
	if ((_disassemble_a > _disassemble_b)); then
		die "ERROR: ${_disassemble_a} (A) is greater than ${_disassemble_b} (B)."
	fi

	local _datasize_array=($(_polap_lib_array-numbers-a2b \
		"$_disassemble_a" \
		"$_disassemble_b" \
		"$_disassemble_n"))
	_polap_log3 "      generated step array (${#_datasize_array[@]}): ${_datasize_array[@]:0:3} ... ${_datasize_array[@]: -2}"

	local _datasize_array_length=${#_datasize_array[@]}

	_polap_log3_cmd rm -f "${_index_table}"
	local i=0
	for ((i = 0; i < _datasize_array_length; i++)); do
		sampling_datasize="${_datasize_array[$i]}"
		printf "%d\t%d\n" $i $sampling_datasize >>"${_index_table}"
	done
}
