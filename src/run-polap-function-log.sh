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
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

# default: _arg_verbose=1
# --quiet: _arg_verbose=0
# -v: _arg_verbose=2
# -v -v: _arg_verbose=3
# -v -v -v: _arg_verbose=4
# Function to handle verbose output
function verbose_echo() {
	local msg_level=$1 # The verbosity level of this message
	shift              # Shift arguments to access the actual message
	local message="$@"

	# Only print if the current verbosity level is greater than or equal to the message level
	if [ "${_arg_verbose}" -ge "$msg_level" ]; then
		# echo "${_arg_verbose} > ${msg_level}: $message"
		echo "$message"
	fi
}

# Function to handle verbose output
function verbose_echo_trim() {
	local msg_level=$1 # The verbosity level of this message
	shift              # Shift arguments to access the actual message
	# local message="$@"

	# Only print if the current verbosity level is greater than or equal to the message level
	if [ "${_arg_verbose}" -ge "$msg_level" ]; then
		# echo "${_arg_verbose} > ${msg_level}: $message"
		local message=$(echo "$@" | tr -d '\n' | tr -s '[:space:]' ' ')
		echo "command: $message"
	fi
}

# Function to handle verbose output
function verbose_echo_no_newline() {
	local msg_level=$1 # The verbosity level of this message
	shift              # Shift arguments to access the actual message
	local message="$@"

	# Only print if the current verbosity level is greater than or equal to the message level
	if [ "${_arg_verbose}" -ge "$msg_level" ]; then
		# echo "${_arg_verbose} > ${msg_level}: $message"
		echo -n "$message"
	fi
}

function verbose_echo_newline() {
	local msg_level=$1 # The verbosity level of this message
	shift              # Shift arguments to access the actual message
	local message="$@"

	# Only print if the current verbosity level is greater than or equal to the message level
	if [ "${_arg_verbose}" -ge "$msg_level" ]; then
		# echo "${_arg_verbose} > ${msg_level}: $message"
		echo ""
		echo "$message"
	fi
}

function verbose_cat() {
	local msg_level=$1 # The verbosity level of this message
	shift              # Shift arguments to access the actual message
	local message="$@"

	# Only print if the current verbosity level is greater than or equal to the message level
	if [ "${_arg_verbose}" -ge "$msg_level" ]; then
		# echo "${_arg_verbose} > ${msg_level}: $message"
		local n=$(wc -l <"${message}")
		echo "FILE ($n lines): $message"
		cat "$message"
	fi
}

function verbose_head() {
	local msg_level=$1 # The verbosity level of this message
	shift              # Shift arguments to access the actual message
	local message="$@"

	# Only print if the current verbosity level is greater than or equal to the message level
	if [ "${_arg_verbose}" -ge "$msg_level" ]; then
		# echo "${_arg_verbose} > ${msg_level}: $message"
		local n=$(wc -l <"${message}")
		echo "FILE head ($n lines): $message"
		head "$message"
		echo "..."
	fi
}

function verbose_column() {
	local msg_level=$1 # The verbosity level of this message
	shift              # Shift arguments to access the actual message
	local message="$@"

	# Only print if the current verbosity level is greater than or equal to the message level
	if [ "${_arg_verbose}" -ge "$msg_level" ]; then
		# echo "${_arg_verbose} > ${msg_level}: $message"
		local n=$(wc -l <"${message}")
		echo "FILE ($n lines): $message"
		column -t "$message"
	fi
}

# default: no stderr output
# print to stderr if --verbose
function echoerr() { verbose_echo 2 "$@" 1>&2; }

# default: output to log but no stderr output
# print to stderr if --verbose
function echoall() {
	verbose_echo 2 "$@" >&3
	verbose_echo 1 "$@"
}

# function yell() { verbose_echo 1 "$0: $*" >&2; }
function yell() {
	verbose_echo 1 "$@" >&2
}

function die() {
	verbose_echo 0 "$@" 1>&2
	verbose_echo 0 "$@"
	exit $EXIT_FAIL
}

function _polap_die() {
	verbose_echo 0 "$@" 1>&2
	verbose_echo 0 "$@"
	exit $EXIT_FAIL
}

# Helper function for logging
# only to the screen with --verbose
function log1_file() {
	echoerr "FILE: $1"
}

function _polap_log0_only() {
	verbose_echo 0 "$@"
}

# --quiet level
# log only to the log file
function _polap_log0() {
	verbose_echo 0 "$@"
	# verbose_echo 1 "$@" 1>&2
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_echo 1 "$@" >&3
	else
		verbose_echo 1 "$@" >&2
	fi
}

# log level 1 to the log file
# log level 0 to the screen
function _polap_log1() {
	verbose_echo 0 "$@"
	# verbose_echo 2 "$@" 1>&2
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_echo 2 "$@" >&3
	else
		verbose_echo 2 "$@" >&2
	fi
}

# log level 2 to the log file
# log level 1 to the screen
function _polap_log2() {
	verbose_echo 0 "$@"
	# verbose_echo 3 "$@" 1>&2
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_echo 3 "$@" >&3
	else
		verbose_echo 3 "$@" >&2
	fi
}

# log level 3 to the log file
# log level 2 to the screen
function _polap_log3() {
	verbose_echo 0 "$@"
	# verbose_echo 4 "$@" 1>&2
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_echo 4 "$@" >&3
	else
		verbose_echo 4 "$@" >&2
	fi
}

function _polap_log0_n() {
	verbose_echo_no_newline 0 "$@"
	verbose_echo_no_newline 1 "$@" >&3
}

function _polap_log0_cmd() {
	verbose_echo_trim 0 "$@"
	verbose_echo_trim 1 "$@" >&3
	"$@"
}

# log level 1 to the log file
# log level 0 to the screen
function _polap_log1_cmd() {
	verbose_echo_trim 0 "$@"
	verbose_echo_trim 2 "$@" >&3
	"$@"
}

# log level 2 to the log file
# log level 1 to the screen
function _polap_log2_cmd() {
	verbose_echo_trim 0 "$@"
	verbose_echo_trim 3 "$@" >&3
	"$@"
}

function _polap_log3_cmd() {
	verbose_echo_trim 0 "$@"
	verbose_echo_trim 4 "$@" >&3
	"$@"
}

function _polap_log2_pipe() {
	verbose_echo_trim 0 "$@"
	verbose_echo_trim 3 "$@" >&3
	eval "$@"
}

function _polap_log3_pipe() {
	verbose_echo_trim 0 "$@"
	verbose_echo_trim 4 "$@" >&3
	eval "$@"
}

function _polap_log0_log() {
	_polap_log0 "LOG: $@"
}

function _polap_log1_log() {
	_polap_log1 "LOG: $@"
}

function _polap_log2_log() {
	_polap_log2 "LOG: $@"
}

function _polap_log3_log() {
	_polap_log3 "LOG: $@"
}

function _polap_echo0() {
	verbose_echo 0 "$@" >&3
}

function _polap_log0_cat() {
	verbose_cat 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_cat 1 "$@" >&3
	else
		verbose_cat 1 "$@" >&2
	fi
}

function _polap_log1_cat() {
	verbose_cat 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_cat 2 "$@" >&3
	else
		verbose_cat 2 "$@" >&2
	fi
}

function _polap_log2_cat() {
	verbose_cat 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_cat 3 "$@" >&3
	else
		verbose_cat 3 "$@" >&2
	fi
}

function _polap_log3_cat() {
	verbose_cat 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_cat 4 "$@" >&3
	else
		verbose_cat 4 "$@" >&2
	fi
}

function _polap_log0_head() {
	verbose_head 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_head 1 "$@" >&3
	else
		verbose_head 1 "$@" >&2
	fi
}

function _polap_log1_head() {
	verbose_head 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_head 2 "$@" >&3
	else
		verbose_head 2 "$@" >&2
	fi
}

function _polap_log2_head() {
	verbose_head 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_head 3 "$@" >&3
	else
		verbose_head 3 "$@" >&2
	fi
}

function _polap_log3_head() {
	verbose_head 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_head 4 "$@" >&3
	else
		verbose_head 4 "$@" >&2
	fi
}

function _polap_log0_column() {
	verbose_column 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_column 1 "$@" >&3
	else
		verbose_column 1 "$@" >&2
	fi
}

function _polap_log1_column() {
	verbose_column 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_column 2 "$@" >&3
	else
		verbose_column 2 "$@" >&2
	fi
}

function _polap_log2_column() {
	verbose_column 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_column 3 "$@" >&3
	else
		verbose_column 3 "$@" >&2
	fi
}

function _polap_log3_column() {
	verbose_column 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_column 4 "$@" >&3
	else
		verbose_column 4 "$@" >&2
	fi
}

function _polap_log0_file() {
	_polap_log0 "FILE: $@"
}

function _polap_log1_file() {
	_polap_log1 "FILE: $@"
}

function _polap_log2_file() {
	_polap_log2 "FILE: $@"
}

function _polap_log3_file() {
	_polap_log3 "FILE: $@"
}

function _polap_log_function() {
	verbose_echo_newline 3 "$@"
	verbose_echo_newline 4 "$@" >&3
}
