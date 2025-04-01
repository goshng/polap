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

function _polap_lib_timing-convert_to_hours_or_minutes {
	local time_string="$1"
	local _hm=0
	local _hours=""

	if [[ "$time_string" =~ ^([0-9]+):([0-9]{2}):([0-9]{2})$ ]]; then
		# h:mm:ss format
		local h="${BASH_REMATCH[1]}"
		local m="${BASH_REMATCH[2]}"
		local s="${BASH_REMATCH[3]}"
		_hm=$(bc <<<"scale=1; $h + $m / 60 + $s / 3600")
		# Use bc for the comparison
		if (($(echo "${_hm} < 1" | bc -l))); then
			_hm=$(bc <<<"scale=1; $m + $s / 60")
			_hours="${_hm}m"
		else
			_hours="${_hm}h"
		fi
	elif [[ "$time_string" =~ ^([0-9]+):([0-9]{2})\.([0-9]{2})$ ]]; then
		# mm:ss.ss format
		local m="${BASH_REMATCH[1]}"
		local s="${BASH_REMATCH[2]}.${BASH_REMATCH[3]}"

		_hm=$(bc <<<"scale=1; $m + $s / 60")
		_hours="${_hm}m"
	else
		echo "Invalid time format: $time_string" >&2
		return 1
	fi

	echo "${_hours}"
}

# read -r memory1 time1 < <(parse_timing Eucalyptus_pauciflora 3)
function _polap_lib_timing-parse-timing {
	local _timing_file="${1}"
	# local _v1="${1}"
	# local j="${2}"
	local _memory_gb
	local _total_hours
	# local _timing_file=${_v1}/timing-${j}.txt

	if [[ -s "${_timing_file}" ]]; then
		local _time_wga=$(grep 'Elapsed' "${_timing_file}" | head -1)
		local _memory_wga=$(grep 'Maximum resident set size' "${_timing_file}" | head -1)
		# Extract the number in kilobytes
		local _memory_kbytes=$(echo "$_memory_wga" | grep -oE "[0-9]+")
		# Convert kilobytes to gigabytes
		local _memory_gb=$(echo "scale=2; $_memory_kbytes / 1048576" | bc)
		# Extract the time portion using grep with regex
		# time_only=$(echo "$_time_wga" | grep -oE "[0-9]+(:[0-9]{2}){1,2}")
		local time_only=$(grep "Elapsed (wall clock) time (h:mm:ss or m:ss):" "$_timing_file" | awk -F': ' '{print $2}')

		local _total_hours=$(_polap_lib_timing-convert_to_hours_or_minutes "${time_only}")
	else
		_memory_gb=0
		_total_hours=0
	fi

	echo "${_memory_gb} ${_total_hours}"
}

# Function to set the start time
_polap_set_start_time() {
	_polap_var_start_time=$(date +%s)
}

# Function to compute elapsed time and return formatted string
_polap_get_elapsed_time() {
	local current_time=$(date +%s)
	local elapsed=$((current_time - _polap_var_start_time))

	if ((elapsed < 60)); then
		echo "${elapsed}s" # Seconds
	elif ((elapsed < 3600)); then
		echo "$((elapsed / 60))m" # Minutes
	elif ((elapsed < 604800)); then
		echo "$((elapsed / 3600))h" # Hours
	else
		echo "$((elapsed / 86400))d" # Days
	fi
}

_polap_get_time_format() {
	local elapsed=$1

	if ((elapsed < 60)); then
		echo "${elapsed}s" # Seconds
	elif ((elapsed < 3600)); then
		echo "$((elapsed / 60))m" # Minutes
	elif ((elapsed < 604800)); then
		echo "$((elapsed / 3600))h" # Hours
	else
		echo "$((elapsed / 86400))d" # Days
	fi
}

function _polap_lib_timing-reset {
	_total_iterations="${1}"
}

# put it at the new start of an iteration
function _polap_lib_timing-set {
	_polap_set_start_time
}

# put it at the end of an iteration
# arg1: the index of a current step, e.g., 0
# arg2: total number of the steps
function _polap_lib_timing-step {
	local i="${1}"
	local _total_iterations="${2}"
	local _weight="${3:-2}"
	local _actual_total_iterations
	local time_per_iteration
	local remaining_iterations
	local remaining_time

	# Progress
	# Calculate elapsed time and remaining time
	time_per_iteration=$(($(date +%s) - _polap_var_start_time))
	if [[ "${_weight}" == "1" ]]; then
		_actual_total_iterations=$(((_total_iterations - i) * (_total_iterations - i + 1) / 2))
	else
		_actual_total_iterations=$((_total_iterations - i))
	fi
	remaining_iterations=$((_actual_total_iterations - 1))
	remaining_time=$((remaining_iterations * time_per_iteration))
	status="    iteration $((i + 1))/${_total_iterations}, remaining time: $(_polap_get_time_format ${remaining_time})"
	echo "$status"
}
