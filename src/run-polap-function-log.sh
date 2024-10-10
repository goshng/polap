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

################################################################################
# for the magic logit function
################################################################################
function logit() {
	while read; do
		# echo "$(date) $REPLY" >> ${LOG_FILE}
		# -1 means "current time"
		# printf "[%(%Y-%m-%d %T)T] %s\n" -1 "$REPLY" >> ${LOG_FILE}
		printf "[%s] %s\n" "$(date +"%Y-%m-%d %T")" "$REPLY" >>${LOG_FILE}
	done
}

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
		head "$message"
	fi
}

# default: no stderr output
# print to stderr if --verbose
function echoerr() { verbose_echo 2 "$@" 1>&2; }

# default: output to log but no stderr output
# print to stderr if --verbose
function echoall() {
	verbose_echo 2 "$@" 1>&2
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

# --quiet level
# log only to the log file
function _polap_log0() {
	verbose_echo 0 "$@"
	verbose_echo 1 "$@" 1>&2
}

# log level 1 to the log file
# log level 0 to the screen
function _polap_log1() {
	verbose_echo 1 "$@"
	verbose_echo 2 "$@" 1>&2
}

# log level 2 to the log file
# log level 1 to the screen
function _polap_log2() {
	verbose_echo 2 "$@"
	verbose_echo 3 "$@" 1>&2
}

# log level 3 to the log file
# log level 2 to the screen
function _polap_log3() {
	verbose_echo 3 "$@"
	verbose_echo 4 "$@" 1>&2
}

function _polap_log3_cmd() {
	verbose_echo 3 "$@"
	verbose_echo 4 "$@" 1>&2
	"$@"
}

function _polap_log3_pipe() {
	verbose_echo 3 "$@"
	verbose_echo 4 "$@" 1>&2
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
	verbose_echo 0 "$@" 1>&2
}

function _polap_log0_cat() {
	verbose_cat 0 "$@"
	verbose_cat 1 "$@" 1>&2
}

function _polap_log1_cat() {
	verbose_cat 1 "$@"
	verbose_cat 2 "$@" 1>&2
}

function _polap_log2_cat() {
	verbose_cat 2 "$@"
	verbose_cat 3 "$@" 1>&2
}

function _polap_log3_cat() {
	verbose_cat 3 "$@"
	verbose_cat 4 "$@" 1>&2
}

function _polap_log_function() {
	verbose_echo_newline 2 "$@"
	verbose_echo_newline 3 "$@" 1>&2
}
