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

function _polap_lib_timing-get_system_info {
	echo "===== System Information ====="

	# Hostname
	echo "Hostname: $(hostname)"

	# Kernel version
	echo "Kernel Version: $(uname -r)"

	# OS version
	if [[ -f /etc/os-release ]]; then
		echo "OS Version: $(source /etc/os-release && echo "$PRETTY_NAME")"
	else
		echo "OS Version: Unknown"
	fi

	# CPU
	echo "CPU Model: $(lscpu | grep 'Model name' | sed 's/Model name:[ \t]*//')"
	echo "CPU Cores: $(nproc)"

	# Memory
	echo "Total Memory: $(free -h | awk '/^Mem:/ {print $2}')"

	# Uptime
	echo "System Uptime: $(uptime -p)"

	# Load average
	echo "Load Average: $(uptime | awk -F'load average: ' '{print $2}')"

	# Disk usage
	echo "Disk Usage (current dir):"
	df -h .

	# Filesystem type
	echo "Filesystem Type (current dir): $(df -T . | tail -1 | awk '{print $2}')"

	# Mount point device
	mount_src=$(df -P . | tail -1 | awk '{print $1}')
	echo "Mount Source: $mount_src"

	# Storage type (HDD/SSD) for current directory
	disk_device=$(lsblk -no pkname "$mount_src" 2>/dev/null)
	if [[ -n "$disk_device" ]]; then
		rotational=$(cat /sys/block/"$disk_device"/queue/rotational 2>/dev/null)
		if [[ "$rotational" == "0" ]]; then
			echo "Storage Type (current dir): SSD"
		elif [[ "$rotational" == "1" ]]; then
			echo "Storage Type (current dir): HDD"
		else
			echo "Storage Type (current dir): Unknown"
		fi
	else
		echo "Storage Device Info: Unknown"
	fi

	# GPU info
	echo "GPU(s):"
	lspci | grep -i 'vga\|3d\|display' || echo "No GPU found or lspci not installed"

	echo "==============================="
}

function _polap_lib_timing-sum_time_values {
	local times=("$@")
	local total_minutes=0

	for t in "${times[@]}"; do
		if [[ "$t" == *h ]]; then
			local num=${t%h}
			local minutes=$(echo "$num * 60" | bc)
		elif [[ "$t" == *m ]]; then
			local num=${t%m}
			local minutes=$num
		fi
		total_minutes=$(echo "$total_minutes + $minutes" | bc)
	done

	if (($(echo "$total_minutes >= 60" | bc -l))); then
		local hours=$(echo "scale=2; $total_minutes / 60" | bc)
		echo "${hours}h"
	else
		echo "${total_minutes}m"
	fi
}

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

# Parse GNU time -v command output.
# read -r memory1 time1 < <(parse_timing Eucalyptus_pauciflora 3)
function _polap_lib_timing-parse-timing {
	local _timing_file="${1}"
	# local _v1="${1}"
	# local j="${2}"
	local _memory_gb
	local _total_hours
	# local _timing_file=${_v1}/timing-${j}.txt

	if [[ -s "${_timing_file}" ]]; then
		local _time_wga=$(grep 'Elapsed' "${_timing_file}" | tail -1)
		local _memory_wga=$(grep 'Maximum resident set size' "${_timing_file}" | tail -1)
		# Extract the number in kilobytes
		local _memory_kbytes=$(echo "$_memory_wga" | grep -oE "[0-9]+")
		# Convert kilobytes to gigabytes
		local _memory_gb=$(echo "scale=2; $_memory_kbytes / 1048576" | bc)
		# Extract the time portion using grep with regex
		# time_only=$(echo "$_time_wga" | grep -oE "[0-9]+(:[0-9]{2}){1,2}")
		local time_only=$(grep "Elapsed (wall clock) time (h:mm:ss or m:ss):" "$_timing_file" | awk -F': ' '{print $2}' | tail -1)

		local _total_hours=$(_polap_lib_timing-convert_to_hours_or_minutes "${time_only}")
	else
		_memory_gb=0
		_total_hours=0
	fi

	echo "${_memory_gb} ${_total_hours}"
}

_polap_lib_timing-parse-multiple-time-only() {
	local _timing_file="$1"
	local max_kb=0
	local count=0

	if [[ -s "$_timing_file" ]]; then
		while IFS= read -r line; do
			local kb=$(echo "$line" | grep -oE '[0-9]+')
			if [[ -n "$kb" ]]; then
				((count++))
				if ((kb > max_kb)); then
					max_kb=$kb
				fi
			fi
		done < <(grep 'Maximum resident set size' "$_timing_file")
	fi

	# Convert to GB with two decimal precision
	local max_gb=$(echo "scale=2; $max_kb / 1048576" | bc)
	echo "$max_gb $count"
}

_polap_lib_timing-parse-cumulative-timing() {
	local _timing_file="$1"
	local total_seconds=0

	if [[ -s "$_timing_file" ]]; then
		while IFS= read -r line; do
			local time_str=$(echo "$line" | awk -F': ' '{print $NF}')
			local seconds=0

			if [[ "$time_str" =~ ^([0-9]+):([0-9]{2}):([0-9]{2})$ ]]; then
				local h="${BASH_REMATCH[1]}"
				local m="${BASH_REMATCH[2]}"
				local s="${BASH_REMATCH[3]}"
				seconds=$((10#$h * 3600 + 10#$m * 60 + 10#$s))
			elif [[ "$time_str" =~ ^([0-9]+):([0-9]{2})\.([0-9]+)$ ]]; then
				local m="${BASH_REMATCH[1]}"
				local s="${BASH_REMATCH[2]}"
				seconds=$((10#$m * 60 + 10#$s))
			else
				continue
			fi

			total_seconds=$((total_seconds + seconds))
		done < <(grep 'Elapsed (wall clock) time' "$_timing_file")

		local total_h=$((total_seconds / 3600))
		local rem=$((total_seconds % 3600))
		local total_m=$((rem / 60))
		local total_s=$((rem % 60))
		local formatted_time=$(printf "%d:%02d:%02d" "$total_h" "$total_m" "$total_s")

		local hm_str=$(_polap_lib_timing-convert_to_hours_or_minutes "$formatted_time")
		echo "$hm_str"
	else
		echo "0h"
	fi
}

# Function to set the start time
_polap_set_start_time() {
	_polap_var_start_time=$(date +%s)
}

# Function to compute elapsed time and return formatted string
_polap_get_elapsed_time() {
	local _brg_start_time="${1:-${_polap_var_start_time}}"
	local current_time=$(date +%s)
	local elapsed=$((current_time - _brg_start_time))

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

function _polap_lib_timing-step-reset {
	_polap_var_total_time=0
}

# put it at the end of an iteration
# arg1: the index of a current step, e.g., 0
# arg2: total number of the steps
function _polap_lib_timing-step {
	local i="${1}"
	local _total_iterations="${2}"
	local _weight="${3:-2}"
	local _stage="${4:-0}"
	local _actual_total_iterations
	local time_per_iteration
	local remaining_iterations
	local remaining_time

	# Progress
	# Calculate elapsed time and remaining time
	time_per_iteration=$(($(date +%s) - _polap_var_start_time))
	# not working in a subshell
	# _polap_var_total_time=$((_polap_var_total_time + time_per_iteration))
	if [[ "${_weight}" == "1" ]]; then
		_actual_total_iterations=$(((_total_iterations - i) * (_total_iterations - i + 1) / 2))
	else
		_actual_total_iterations=$((_total_iterations - i))
	fi
	remaining_iterations=$((_actual_total_iterations - 1))
	remaining_time=$((remaining_iterations * time_per_iteration))
	status="  ${_stage} iteration $((i + 1))/${_total_iterations}, remaining time: $(_polap_get_time_format ${remaining_time})"
	echo "$status"
}
