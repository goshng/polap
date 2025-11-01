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
# Polap writes activities in <out>/polap.log file.
# This script defines various log functions.
# See Also:
# run-polap-function-template.sh
# TODO: rename: polap-lib-log.sh
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

# default: _arg_verbose=1
# --quiet: _arg_verbose=0
# -v: _arg_verbose=2
# -v -v: _arg_verbose=3
# -v -v -v: _arg_verbose=4
# Function to handle verbose output
function verbose_echo {
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
function verbose_echo_trim {
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
function verbose_echo_no_newline {
	local msg_level=$1 # The verbosity level of this message
	shift              # Shift arguments to access the actual message
	local message="$@"

	# Only print if the current verbosity level is greater than or equal to the message level
	if [ "${_arg_verbose}" -ge "$msg_level" ]; then
		# echo "${_arg_verbose} > ${msg_level}: $message"
		echo -n "$message"
	fi
}

# Function to handle verbose output
function verbose_echo_no_newline_ne {
	local msg_level=$1 # The verbosity level of this message
	shift              # Shift arguments to access the actual message
	local message="$@"

	# Only print if the current verbosity level is greater than or equal to the message level
	if [ "${_arg_verbose}" -eq "$msg_level" ]; then
		# echo "${_arg_verbose} > ${msg_level}: $message"
		echo -ne "$message"
	fi
}

function verbose_echo_newline {
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

function verbose_cat {
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

function verbose_head {
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

function verbose_column {
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
function echoerr { verbose_echo 2 "$@" 1>&2; }

# default: output to log but no stderr output
# print to stderr if --verbose
function echoall {
	verbose_echo 2 "$@" >&3
	verbose_echo 1 "$@"
}

# function yell { verbose_echo 1 "$0: $*" >&2; }
function yell {
	verbose_echo 1 "$@" >&2
}

function die {
	verbose_echo 0 "$@" 1>&2
	verbose_echo 0 "$@"
	exit $EXIT_FAIL
}

function _polap_die {
	verbose_echo 0 "$@" 1>&2
	verbose_echo 0 "$@"
	exit $EXIT_FAIL
}

# Helper function for logging
# only to the screen with --verbose
function log1_file {
	echoerr "FILE: $1"
}

function _polap_log0_only {
	verbose_echo 0 "$@"
}

# --quiet level
# log only to the log file
#
# if [[ "${_POLAP_RELEASE}" == "1" ]]; then
# 	function _polap_log0 {
# 		verbose_echo 0 "$@"
# 		# verbose_echo 1 "$@" 1>&2
# 		if [[ "${_arg_log_stderr}" = "off" ]]; then
# 			verbose_echo 1 "$@" >&3
# 		else
# 			verbose_echo 1 "$@" >&2
# 		fi
# 	}
# else
# 	function _polap_log0 {
# 		local func="${FUNCNAME[1]}"
# 		local file="$(basename ${BASH_SOURCE[1]})"
# 		local line="${BASH_LINENO[0]}"
# 		local tag="[$func@$file:$line]"

# 		verbose_echo 0 "${tag} $@"

# 		if [[ "${_arg_log_stderr}" = "off" ]]; then
# 			verbose_echo 1 "${tag} $@" >&3
# 		else
# 			verbose_echo 1 "${tag} $@" >&2
# 		fi
# 	}
# fi

# log level 1 to the log file
# log level 0 to the screen
# function _polap_log1 {
# 	verbose_echo 0 "$@"
# 	# verbose_echo 2 "$@" 1>&2
# 	if [[ "${_arg_log_stderr}" = "off" ]]; then
# 		verbose_echo 2 "$@" >&3
# 	else
# 		verbose_echo 2 "$@" >&2
# 	fi
# }

# log level 2 to the log file
# log level 1 to the screen
# function _polap_log2 {
# 	verbose_echo 0 "$@"
# 	# verbose_echo 3 "$@" 1>&2
# 	if [[ "${_arg_log_stderr}" = "off" ]]; then
# 		verbose_echo 3 "$@" >&3
# 	else
# 		verbose_echo 3 "$@" >&2
# 	fi
# }

# log level 3 to the log file
# log level 2 to the screen
# function _polap_log3 {
# 	verbose_echo 0 "$@"
# 	# verbose_echo 4 "$@" 1>&2
# 	if [[ "${_arg_log_stderr}" = "off" ]]; then
# 		verbose_echo 4 "$@" >&3
# 	else
# 		verbose_echo 4 "$@" >&2
# 	fi
# }

# --quiet level
# log only to the log file

if [[ "${_POLAP_RELEASE}" == "1" ]]; then
	function _polap_log0 {
		verbose_echo 0 "$@"
		if [[ "${_arg_log_stderr}" = "off" ]]; then
			verbose_echo 1 "$@" >&3
		else
			verbose_echo 1 "$@" >&2
		fi
	}
	function _polap_log1 {
		verbose_echo 0 "$@"
		if [[ "${_arg_log_stderr}" = "off" ]]; then
			verbose_echo 2 "$@" >&3
		else
			verbose_echo 2 "$@" >&2
		fi
	}
	function _polap_log2 {
		verbose_echo 0 "$@"
		if [[ "${_arg_log_stderr}" = "off" ]]; then
			verbose_echo 3 "$@" >&3
		else
			verbose_echo 3 "$@" >&2
		fi
	}
	function _polap_log3 {
		verbose_echo 0 "$@"
		if [[ "${_arg_log_stderr}" = "off" ]]; then
			verbose_echo 4 "$@" >&3
		else
			verbose_echo 4 "$@" >&2
		fi
	}
else
	# Debug mode with func@file:line tags
	__polap_log_with_tag() {
		local level="$1"
		shift
		local outlevel="$1"
		shift
		local func="${FUNCNAME[2]}"
		local file="$(basename "${BASH_SOURCE[2]}")"
		local line="${BASH_LINENO[1]}"
		local tag="[$func@$file:$line]"

		verbose_echo 0 "${tag} $*"
		if [[ "${_arg_log_stderr}" = "off" ]]; then
			verbose_echo "${outlevel}" "${tag} $*" >&3
		else
			verbose_echo "${outlevel}" "${tag} $*" >&2
		fi
	}

	function _polap_log0 { __polap_log_with_tag 0 1 "$@"; }
	function _polap_log1 { __polap_log_with_tag 0 2 "$@"; }
	function _polap_log2 { __polap_log_with_tag 0 3 "$@"; }
	function _polap_log3 { __polap_log_with_tag 0 4 "$@"; }
fi

if [[ "${_POLAP_DEBUG}" == "1" ]]; then
	function _polap_log0_dev {
		local func="${FUNCNAME[1]}"
		local file="$(basename ${BASH_SOURCE[1]})"
		local line="${BASH_LINENO[0]}"
		local tag="[$func@$file:$line]"

		verbose_echo 0 "${tag} $@"

		if [[ "${_arg_log_stderr}" = "off" ]]; then
			verbose_echo 1 "${tag} $@" >&3
		else
			verbose_echo 1 "${tag} $@" >&2
		fi
	}
else
	function _polap_log0_dev {
		:
	}
fi

function _polap_log0_n {
	verbose_echo_no_newline 0 "$@"
	verbose_echo_no_newline 1 "$@" >&3
}

function _polap_log0_ne {
	# verbose_echo_no_newline_ne 0 "$@"
	verbose_echo_no_newline_ne 1 "$@" >&3
}

function _polap_log1_ne {
	# verbose_echo_no_newline_ne 0 "$@"
	verbose_echo_no_newline_ne 2 "$@" >&3
}

# function _polap_log0_cmd {
# 	verbose_echo_trim 0 "$@"
# 	verbose_echo_trim 1 "$@" >&3
# 	"$@"
# }
#
# # log level 1 to the log file
# # log level 0 to the screen
# function _polap_log1_cmd {
# 	verbose_echo_trim 0 "$@"
# 	verbose_echo_trim 2 "$@" >&3
# 	"$@"
# }
#
# # log level 2 to the log file
# # log level 1 to the screen
# function _polap_log2_cmd {
# 	verbose_echo_trim 0 "$@"
# 	verbose_echo_trim 3 "$@" >&3
# 	"$@"
# }
#
# function _polap_log3_cmd {
# 	verbose_echo_trim 0 "$@"
# 	verbose_echo_trim 4 "$@" >&3
# 	"$@"
# }

# Base runner: logs to screen at level 0 and to FD3 at the level you pass.
_polap_log_cmd() {
	local _fd3_level="$1"
	shift
	verbose_echo_trim 0 "$@" || true
	verbose_echo_trim "${_fd3_level}" "$@" >&3 || true
	"$@"
}

# Thin wrappers (only the FD3 verbosity differs)
_polap_log0_cmd() { _polap_log_cmd 1 "$@"; }
_polap_log1_cmd() { _polap_log_cmd 2 "$@"; }
_polap_log2_cmd() { _polap_log_cmd 3 "$@"; }
_polap_log3_cmd() { _polap_log_cmd 4 "$@"; }

# function _polap_log0_pipe {
# 	verbose_echo_trim 0 "$@"
# 	verbose_echo_trim 1 "$@" >&3
# 	eval "$@"
# }
#
# function _polap_log1_pipe {
# 	verbose_echo_trim 0 "$@"
# 	verbose_echo_trim 2 "$@" >&3
# 	eval "$@"
# }
#
# function _polap_log2_pipe {
# 	verbose_echo_trim 0 "$@"
# 	verbose_echo_trim 3 "$@" >&3
# 	eval "$@"
# }
#
# function _polap_log3_pipe {
# 	verbose_echo_trim 0 "$@"
# 	verbose_echo_trim 4 "$@" >&3
# 	eval "$@"
# }

#!/usr/bin/env bash
# Version: v0.1.0

# Base: log to screen at level 0, to FD3 at the level you pass.
# _soft=0 → normal (return real status)
# _soft=1 → "no-exit" (suspend -e & ERR trap; always return 0, stash raw in _POLAP_LOG3_LAST_STATUS)
_polap_log_pipe_impl() {
	local _fd3_level="$1" _soft="$2"
	shift 2
	local cmd="$*"

	verbose_echo_trim 0 "$cmd" || true
	verbose_echo_trim "${_fd3_level}" "$cmd" >&3 || true

	if ((_soft)); then
		# --- soft mode: don't trip errexit even on non-zero ---
		local had_e=0 st saved_errtrap=""
		case $- in *e*) had_e=1 ;; esac
		saved_errtrap="$(trap -p ERR || true)"
		set +e
		trap - ERR
		eval "$cmd"
		st=$?
		[[ -n "$saved_errtrap" ]] && eval "$saved_errtrap"
		((had_e)) && set -e
		_POLAP_LOG3_LAST_STATUS=$st
		# (optional) warn on non-zero:
		((st != 0)) && printf '[%(%F %T)T WARN] non-zero exit %d: %s\n' -1 "$st" "$cmd" >&2
		return 0
	else
		# --- normal mode: return real status (respects set -e at call site) ---
		eval "$cmd"
	fi
}

# Thin wrappers (only FD3 verbosity and softness differ)
_polap_log0_pipe() { _polap_log_pipe_impl 1 0 "$@"; }
_polap_log1_pipe() { _polap_log_pipe_impl 2 0 "$@"; }
_polap_log2_pipe() { _polap_log_pipe_impl 3 0 "$@"; }
_polap_log3_pipe() { _polap_log_pipe_impl 4 0 "$@"; } # Avoid err-exit

_polap_log3_pipe_avoid_errexit() { _polap_log_pipe_impl 4 1 "$@"; }

# Version: v0.4.1
# Features:
#   - --argv: run without eval (safer)
#   - --ok "0 1": treat listed exit codes as success (e.g., grep 1 = no match)
#   - --out FILE: redirect stdout to FILE (truncate)
#   - --append FILE: redirect stdout to FILE (append)
#   - --status-var VAR: store raw status (before ok-coercion) in VAR
#   - logs the redirection text in the display line, too.
#
# _polap_log_pipe_no_exit() — safe runner with logging & selective OK exit-codes.
# Wrappers _polap_log{0,1,2,3}_pipe_no_exit() only change the FD3 log level.
#
_polap_log_pipe_no_exit() {
	local ok_codes="" mode="string" status_var="" out_path="" append_path=""
	local __fd3_level=4 # default for “log3”; wrappers override via --_fd3-level
	while [[ $# -gt 0 ]]; do
		case "$1" in
		--ok | --ok-codes)
			ok_codes="$2"
			shift 2
			;;
		--status-var)
			status_var="$2"
			shift 2
			;;
		--argv)
			mode="argv"
			shift
			;;
		--out)
			out_path="$2"
			shift 2
			;;
		--append)
			append_path="$2"
			shift 2
			;;
		--_fd3-level)
			__fd3_level="$2"
			shift 2
			;; # internal
		--)
			shift
			break
			;;
		*) break ;;
		esac
	done

	# Build display string (quoted argv + redirection hint)
	local display s=""
	if [[ "$mode" == "argv" ]]; then
		for a in "$@"; do printf -v s '%s %q' "$s" "$a"; done
		display="${s# }"
		[[ -n "$out_path" ]] && display+=" > ${out_path}"
		[[ -n "$append_path" ]] && display+=" >> ${append_path}"
	else
		display="$1"
	fi

	# Primary line (level 0) and FD3 line at chosen level
	verbose_echo_trim 0 "$display"
	verbose_echo_trim "${__fd3_level}" "$display" >&3 || true

	# Remember state and suspend -e + ERR trap during execution
	local had_e=0 st saved_errtrap=""
	case $- in *e*) had_e=1 ;; esac
	saved_errtrap="$(trap -p ERR || true)"
	set +e
	trap - ERR

	if [[ "$mode" == "argv" ]]; then
		# Ensure parent dir if redirecting
		if [[ -n "$out_path" ]]; then mkdir -p -- "$(dirname -- "$out_path")"; fi
		if [[ -n "$append_path" ]]; then mkdir -p -- "$(dirname -- "$append_path")"; fi

		if [[ -n "$out_path" ]]; then
			"$@" >"$out_path"
			st=$?
		elif [[ -n "$append_path" ]]; then
			"$@" >>"$append_path"
			st=$?
		else
			"$@"
			st=$?
		fi
	else
		eval "$1"
		st=$?
	fi

	# Restore trap and errexit
	[[ -n "$saved_errtrap" ]] && eval "$saved_errtrap"
	((had_e)) && set -e

	# Report raw status if requested
	[[ -n "$status_var" ]] && printf -v "$status_var" '%s' "$st"

	# Whitelist OK codes (e.g., grep: 0 = match, 1 = no match)
	if [[ -n "$ok_codes" ]]; then
		for ok in $ok_codes; do
			[[ "$st" -eq "$ok" ]] && return 0
		done
	fi
	return "$st"
}

# Thin wrappers (only FD3 verbosity differs)
_polap_log0_pipe_no_exit() { _polap_log_pipe_no_exit --_fd3-level 1 "$@"; }
_polap_log1_pipe_no_exit() { _polap_log_pipe_no_exit --_fd3-level 2 "$@"; }
_polap_log2_pipe_no_exit() { _polap_log_pipe_no_exit --_fd3-level 3 "$@"; }
_polap_log3_pipe_no_exit() { _polap_log_pipe_no_exit --_fd3-level 4 "$@"; }

v1_polap_log3_pipe_no_exit() {
	local ok_codes="" mode="string" status_var="" out_path="" append_path=""
	while [[ $# -gt 0 ]]; do
		case "$1" in
		--ok | --ok-codes)
			ok_codes="$2"
			shift 2
			;;
		--status-var)
			status_var="$2"
			shift 2
			;;
		--argv)
			mode="argv"
			shift
			;;
		--out)
			out_path="$2"
			shift 2
			;; # '>'
		--append)
			append_path="$2"
			shift 2
			;; # '>>'
		--)
			shift
			break
			;;
		*) break ;;
		esac
	done

	# Build display string (quoted argv + redir hint)
	local display s=""
	if [[ "$mode" == "argv" ]]; then
		for a in "$@"; do printf -v s '%s %q' "$s" "$a"; done
		display="${s# }"
		if [[ -n "$out_path" ]]; then display+=" > ${out_path}"; fi
		if [[ -n "$append_path" ]]; then display+=" >> ${append_path}"; fi
	else
		display="$1"
	fi

	verbose_echo_trim 0 "$display"
	verbose_echo_trim 4 "$display" >&3 || true

	# Remember -e and current ERR trap; suspend both while running
	local had_e=0 st saved_errtrap=""
	case $- in *e*) had_e=1 ;; esac
	saved_errtrap="$(trap -p ERR || true)"
	set +e
	trap - ERR

	if [[ "$mode" == "argv" ]]; then
		# ensure parent dir if redirecting
		if [[ -n "$out_path" ]]; then install -d -- "$(dirname -- "$out_path")"; fi
		if [[ -n "$append_path" ]]; then install -d -- "$(dirname -- "$append_path")"; fi

		if [[ -n "$out_path" ]]; then
			"$@" >"$out_path"
			st=$?
		elif [[ -n "$append_path" ]]; then
			"$@" >>"$append_path"
			st=$?
		else
			"$@"
			st=$?
		fi
	else
		eval "$1"
		st=$?
	fi

	# restore state
	[[ -n "$saved_errtrap" ]] && eval "$saved_errtrap"
	((had_e)) && set -e

	# raw status out
	[[ -n "$status_var" ]] && printf -v "$status_var" '%s' "$st"

	# whitelist OK codes (e.g., grep: 0=match, 1=no match)
	if [[ -n "$ok_codes" ]]; then
		for ok in $ok_codes; do
			[[ "$st" -eq "$ok" ]] && return 0
		done
	fi
	return "$st"
}

function _polap_log3_pipe_command {
	local cmd

	if [[ -t 0 ]]; then
		# stdin is a terminal -> called with arguments like: _polap_log3_pipe_command ls -l
		cmd="$*"
	else
		# stdin is a pipe -> called with heredoc: _polap_log3_pipe_command <<'EOF' ... EOF
		cmd="$(cat)"
	fi

	verbose_echo_trim 0 "$cmd"
	verbose_echo_trim 4 "$cmd" >&3

	eval "$cmd"
}

# function _polap_log3_cmdout {
# 	local cmd=("$@")
# 	local full_cmd="${cmd[*]} >${_polap_output_dest} 2>&1"
#
# 	verbose_echo_trim 0 "$full_cmd"
# 	verbose_echo_trim 4 "$full_cmd" >&3
#
# 	"${cmd[@]}" >"${_polap_output_dest}" 2>&1
# }

# Run a command and log it at level N (0..3), mirroring _polap_logN behavior.
# Usage:
#   _polap_log_cmdout N cmd arg1 arg2 ...
#   (or use wrappers: _polap_log0_cmdout ... _polap_log3_cmdout)

function _polap_log_cmdout {
	local level="$1"
	shift
	local cmd=("$@")
	local full_cmd="${cmd[*]} >${_polap_output_dest} 2>&1"
	local outlevel=$((level + 1))

	if [[ "${_POLAP_RELEASE}" == "1" ]]; then
		# Release mode: plain
		verbose_echo_trim 0 "$full_cmd"
		if [[ "${_arg_log_stderr}" = "off" ]]; then
			verbose_echo_trim "${outlevel}" "$full_cmd" >&3
		else
			verbose_echo_trim "${outlevel}" "$full_cmd" >&2
		fi
	else
		# Debug mode: tag with func@file:line
		local caller_func="${FUNCNAME[1]}"
		local caller_file="$(basename "${BASH_SOURCE[1]}")"
		local caller_line="${BASH_LINENO[0]}"
		local tag="[$caller_func@${caller_file}:${caller_line}]"

		verbose_echo_trim 0 "${tag} $full_cmd"
		if [[ "${_arg_log_stderr}" = "off" ]]; then
			verbose_echo_trim "${outlevel}" "${tag} $full_cmd" >&3
		else
			verbose_echo_trim "${outlevel}" "${tag} $full_cmd" >&2
		fi
	fi

	# Execute, send stdout+stderr to destination
	"${cmd[@]}" >"${_polap_output_dest}" 2>&1
}

# Compatibility wrappers
_polap_log0_cmdout() { _polap_log_cmdout 0 "$@"; }
_polap_log1_cmdout() { _polap_log_cmdout 1 "$@"; }
_polap_log2_cmdout() { _polap_log_cmdout 2 "$@"; }
_polap_log3_cmdout() { _polap_log_cmdout 3 "$@"; }

function _polap_log0_log {
	_polap_log0 "LOG: $@"
}

function _polap_log1_log {
	_polap_log1 "LOG: $@"
}

function _polap_log2_log {
	_polap_log2 "LOG: $@"
}

function _polap_log3_log {
	_polap_log3 "LOG: $@"
}

function _polap_echo0 {
	verbose_echo 0 "$@" >&3
}

function _polap_log0_cat {
	verbose_cat 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_cat 1 "$@" >&3
	else
		verbose_cat 1 "$@" >&2
	fi
}

function _polap_log1_cat {
	verbose_cat 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_cat 2 "$@" >&3
	else
		verbose_cat 2 "$@" >&2
	fi
}

function _polap_log2_cat {
	verbose_cat 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_cat 3 "$@" >&3
	else
		verbose_cat 3 "$@" >&2
	fi
}

function _polap_log3_cat {
	verbose_cat 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_cat 4 "$@" >&3
	else
		verbose_cat 4 "$@" >&2
	fi
}

function _polap_log0_head {
	verbose_head 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_head 1 "$@" >&3
	else
		verbose_head 1 "$@" >&2
	fi
}

function _polap_log1_head {
	verbose_head 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_head 2 "$@" >&3
	else
		verbose_head 2 "$@" >&2
	fi
}

function _polap_log2_head {
	verbose_head 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_head 3 "$@" >&3
	else
		verbose_head 3 "$@" >&2
	fi
}

function _polap_log3_head {
	verbose_head 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_head 4 "$@" >&3
	else
		verbose_head 4 "$@" >&2
	fi
}

function _polap_log0_column {
	if [[ -s "$1" ]]; then
		verbose_column 0 "$@"
		if [[ "${_arg_log_stderr}" = "off" ]]; then
			verbose_column 1 "$@" >&3
		else
			verbose_column 1 "$@" >&2
		fi
	else
		_polap_log0 "ERROR: no such file: $1"
	fi
}

function _polap_log1_column {
	verbose_column 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_column 2 "$@" >&3
	else
		verbose_column 2 "$@" >&2
	fi
}

function _polap_log2_column {
	verbose_column 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_column 3 "$@" >&3
	else
		verbose_column 3 "$@" >&2
	fi
}

function _polap_log3_column {
	verbose_column 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_column 4 "$@" >&3
	else
		verbose_column 4 "$@" >&2
	fi
}

function _polap_log0_file {
	_polap_log0 "FILE: $@"
}

function _polap_log1_file {
	_polap_log1 "FILE: $@"
}

function _polap_log2_file {
	_polap_log2 "FILE: $@"
}

function _polap_log3_file {
	_polap_log3 "FILE: $@"
}

function _polap_log_function {
	verbose_echo_newline 0 "$@"
	if [[ "${_arg_log_stderr}" = "off" ]]; then
		verbose_echo_newline 4 "$@" >&3
	else
		verbose_echo_newline 4 "$@" >&2
	fi
}

_polap_log_file_exist() {
	local level="$1"
	local msg="$2"
	local file="$3"

	if [[ -s "$file" ]]; then
		local msgout="$msg: $file"
		_polap_log"${level}" "$msgout"
	else
		local msgout="[WARN] no such file exists: $msg: $file"
		_polap_log"${level}" "$msgout"
	fi

}

_polap_log3_file_exist() {
	_polap_log_file_exist 3 "$@"
}

_polap_log2_file_exist() {
	_polap_log_file_exist 2 "$@"
}
_polap_log1_file_exist() {
	_polap_log_file_exist 1 "$@"
}

_polap_log0_file_exist() {
	_polap_log_file_exist 0 "$@"
}
