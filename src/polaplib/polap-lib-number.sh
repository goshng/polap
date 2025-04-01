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

# Function to round a floating-point number
function _polap_lib_number-round-number {
	local number="$1"
	# Add 0.5 before truncating to simulate rounding
	echo "$number + 0.5" | bc | cut -d'.' -f1
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
_polap_lib_number-minimum-two-values-in-files() {
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
