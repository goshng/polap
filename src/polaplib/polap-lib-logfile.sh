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
# Convert numbers between different units.
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

# Print timestamp to both terminal and log file if enabled
_polap_lib_logfile-clock() {
	local clock_flag="${_arg_clock:-off}"
	local log_file="${LOG_FILE:-}"
	local timestamp

	# Only proceed if _arg_clock is 'on'
	[[ "$clock_flag" != "on" ]] && return 0

	timestamp="$(date +"%Y-%m-%d %H:%M:%S")"

	# If LOG_FILE is defined
	if [[ -n "$log_file" ]]; then
		# Case 1: File exists and is writable
		if [[ -w "$log_file" ]]; then
			echo "$timestamp" | tee -a "$log_file" >&1
			return 0

		# Case 2: File doesn't exist but directory is writable
		elif [[ ! -e "$log_file" && -w "$(dirname "$log_file")" ]]; then
			echo "$timestamp" | tee -a "$log_file" >&1
			return 0

		# Case 3: File or directory not writable → fallback to screen
		else
			echo "$timestamp"
			return 0
		fi
	fi

	# If LOG_FILE unset → print only to screen
	echo "$timestamp"
}

# Write POLAP version and executed command to log file (and screen fallback)
_polap_lib_logfile-cmd() {
	local cmd="$1"
	local begin="${2:-BEGIN}"
	local log_file="${LOG_FILE:-}"
	local version="${_polap_version:-unknown}"

	# If LOG_FILE is defined and writable (or creatable), log there
	if [[ -n "$log_file" ]]; then
		if [[ -w "$log_file" || (! -e "$log_file" && -w "$(dirname "$log_file")") ]]; then
			{
				if [[ "${begin}" == "BEGIN" ]]; then
					echo "POLAP: ${version}"
				fi
				echo "CMD: ${cmd}"
			} >>"$log_file"
			return 0
		fi
	fi
}

# Log all variables beginning with _arg_ to the log file (fallback to screen if needed)
_polap_lib_logfile-args() {
	local log_file="${LOG_FILE:-}"

	# Disable 'set -u' temporarily to avoid errors on unset vars
	set +u
	{
		for var in $(compgen -v _arg_); do
			printf "%s=%s\n" "$var" "${!var}"
		done
	} | {
		# Write to log file if valid, otherwise print to screen
		if [[ -n "$log_file" ]]; then
			if [[ -w "$log_file" || (! -e "$log_file" && -w "$(dirname "$log_file")") ]]; then
				cat >>"$log_file"
				set -u
				return 0
			fi
		fi
		cat
	}
	set -u
}

# Log elapsed runtime (in hrs/min/sec) with hostname and command context
_polap_lib_logfile-elapsed() {
	local clock_flag="${_arg_clock:-off}"
	local log_file="${LOG_FILE:-}"
	local elapsed

	# Build elapsed time string using Bash SECONDS
	elapsed="Time at $(hostname): $((SECONDS / 3600))hrs $(((SECONDS / 60) % 60))min $((SECONDS % 60))sec"

	if [[ "$clock_flag" == "on" ]]; then
		# Otherwise, print only to screen
		echo "$elapsed"
	fi

	# If LOG_FILE is valid, append and also print to screen
	if [[ -n "$log_file" ]]; then
		if [[ -w "$log_file" || (! -e "$log_file" && -w "$(dirname "$log_file")") ]]; then
			echo "$elapsed" >>"$log_file"
			return 0
		fi
	fi

}
