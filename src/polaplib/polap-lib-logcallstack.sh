#!/usr/bin/env bash
# polap-lib-logcallstack.sh
# Version: v0.7.0
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
#
# Purpose
#   Unified DEBUG/ERR trapping with a small in-memory ring buffer. When a
#   command fails (non-zero exit with `set -e` + `errtrace`), we print:
#     - FAILED CMD: the exact command (prefer the preview saved by the runner)
#     - the *issuing* site (file:line) and a short “called from” chain
#     - a bash `caller` traceback for extra context
#
# Design
#   - FD3 must be opened by the caller (polap.sh): exec 3>>"$LOG_FILE"
#   - Screen output is gated by _arg_verbose (0..4) and sent to stdout
#     unless _arg_log_stderr=on (then to stderr).
#   - We don’t depend on polap-lib-log.sh internals; this file is
#     self-contained and safe to source early.
#
################################################################################

################################################################################
# Source-once guard
################################################################################
# shellcheck disable=SC1091
if [[ -z "${_POLAPLIB_DIR:-}" ]]; then
	echo "[ERROR] _POLAPLIB_DIR is not set before sourcing polap-lib-logcallstack.sh" >&2
	return 1 2>/dev/null || exit 1
fi
source "${_POLAPLIB_DIR}/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u
if [[ -n "${!_POLAP_INCLUDE_}" ]]; then
	set -u
	return 0
fi
set -u
declare "$_POLAP_INCLUDE_=1"

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
	echo "[ERROR] This script must be sourced, not executed: use 'source $BASH_SOURCE'" >&2
	return 1 2>/dev/null || exit 1
fi

: "${_POLAP_DEBUG:=0}"
: "${_POLAP_RELEASE:=0}"
: "${_arg_verbose:=1}"       # 0=quiet, 1=default, 2=info, 3=detail, 4=debug
: "${_arg_log_stderr:=off}"  # 'on' → screen messages go to STDERR
: "${_POLAP_DBG_MAX:=32}"    # ring size
: "${_POLAP_MAX_CALLERS:=8}" # how many wrapper frames to show
: "${_POLAP_WRAPPER_REGEX:=^(bash|sh|python3?|Rscript)([[:space:]]|$)}"

# Optional xtrace if _POLAP_DEBUG=1 (kept minimal)
if [[ "${_POLAP_DEBUG}" -eq 1 ]]; then
	export PS4='[${EPOCHREALTIME} ${FUNCNAME[0]}@${BASH_SOURCE##*/}:${LINENO}] '
fi

################################################################################
# Small emission helpers (FD3 + screen gating)
################################################################################
_polap__screen_fd() { [[ "${_arg_log_stderr}" == "on" ]] && echo 2 || echo 1; }
_polap__screen_allowed() {
	local lvl="${1:-0}"
	((_arg_verbose >= lvl + 1))
}
_polap__now() { date "+%Y-%m-%d %H:%M:%S"; }

# Print one line to FD3 (always) and to screen if allowed under verbosity
# Args: level  message...
_polap__emit() {
	local level="${1:-0}"
	shift || true
	local line="$*"
	printf "%s\n" "${line}" >&3 || true
	if _polap__screen_allowed "${level}"; then
		local fd
		fd="$(_polap__screen_fd)"
		printf "%s\n" "${line}" >&"$fd" || true
	fi
}

################################################################################
# DEBUG ring buffer
################################################################################
declare -a _POLAP_DBG_CMD=() _POLAP_DBG_FILE=() _POLAP_DBG_LINE=() _POLAP_DBG_FUNC=()
_POLAP_DBG_IDX=0

# Capture the *next* command about to execute
_polap_on_debug() {
	# The simple command to be executed:
	local cmd="${BASH_COMMAND-}"
	# caller 0 → "<line> <func> <file>"
	local __line __func __file
	if read -r __line __func __file < <(caller 0); then
		: "${__func:=main}"
	else
		__file="${BASH_SOURCE[1]-${BASH_SOURCE[0]}}"
		__line="${BASH_LINENO[0]-${LINENO-0}}"
		__func="${FUNCNAME[1]-main}"
	fi

	local idx="${_POLAP_DBG_IDX}"
	_POLAP_DBG_CMD[$idx]="$cmd"
	_POLAP_DBG_FILE[$idx]="${__file}"
	_POLAP_DBG_LINE[$idx]="${__line}"
	_POLAP_DBG_FUNC[$idx]="${__func}"
	_POLAP_DBG_IDX=$(((idx + 1) % _POLAP_DBG_MAX))
}

# Decide whether a command is a “wrapper”
_polap_is_wrapper_cmd() {
	[[ "$1" =~ ${_POLAP_WRAPPER_REGEX} ]]
}

# Walk ring backwards to find most relevant (non-wrapper) command
# Start at head = last written - 1
_polap_ring_find_last_relevant() {
	local head="$(((_POLAP_DBG_IDX + _POLAP_DBG_MAX - 1) % _POLAP_DBG_MAX))"
	local steps=0 i="$head"
	while ((steps < _POLAP_DBG_MAX)); do
		local c="${_POLAP_DBG_CMD[$i]-}"
		if [[ -z "$c" ]]; then
			i=$(((i + _POLAP_DBG_MAX - 1) % _POLAP_DBG_MAX))
			((steps++))
			continue
		fi
		if _polap_is_wrapper_cmd "$c"; then
			i=$(((i + _POLAP_DBG_MAX - 1) % _POLAP_DBG_MAX))
			((steps++))
			continue
		fi
		printf "%s" "$i"
		return 0
	done
	# fallback: return head
	printf "%s" "$head"
}

################################################################################
# Pretty printers
################################################################################
_polap_print_failed_cmd_at() {
	local idx="$1"
	local cmd="${_POLAP_DBG_CMD[$idx]-}"
	local file="${_POLAP_DBG_FILE[$idx]-}"
	local line="${_POLAP_DBG_LINE[$idx]-0}"
	local func="${_POLAP_DBG_FUNC[$idx]-main}"

	local disp=""
	if [[ -n "${POLAP_LAST_CMD_PREVIEW:-}" ]]; then
		disp="${POLAP_LAST_CMD_PREVIEW}"
	elif [[ -n "$cmd" ]]; then
		disp="$cmd"
	else
		disp="${BASH_COMMAND:-}"
	fi

	local ts
	ts="$(_polap__now)"
	local base="${file##*/}"
	# FAILED CMD line
	_polap__emit 0 "[$ts ${func}@${base}:${line}] FAILED CMD: ${disp//$'\n'/\\n}"

	# If we know the *issuing* site (from the runner), print it once
	if [[ -n "${POLAP_LAST_CMD_FILE:-}" && -n "${POLAP_LAST_CMD_LINE:-}" ]]; then
		_polap__emit 0 "[$ts wrapper@${POLAP_LAST_CMD_FILE}:${POLAP_LAST_CMD_LINE}] ISSUED: ${POLAP_LAST_CMD_PREVIEW:-}"
	fi
}

_polap_print_called_from_chain() {
	local frozen_head="$1"
	local ts
	ts="$(_polap__now)"
	local j="$frozen_head" printed=0

	while ((printed < _POLAP_MAX_CALLERS)); do
		local wcmd="${_POLAP_DBG_CMD[$j]-}"
		[[ -z "$wcmd" ]] && break
		if _polap_is_wrapper_cmd "$wcmd"; then
			local wfile="${_POLAP_DBG_FILE[$j]-}"
			local wline="${_POLAP_DBG_LINE[$j]-0}"
			local wfunc="${_POLAP_DBG_FUNC[$j]-main}"
			_polap__emit 0 "[$ts ${wfunc}@${wfile##*/}:${wline}] CALLED FROM: ${wcmd//$'\n'/\\n}"
			printed=$((printed + 1))
			j=$(((j + _POLAP_DBG_MAX - 1) % _POLAP_DBG_MAX))
		else
			break
		fi
	done
}

_polap_print_bash_caller_trace() {
	# Print a mini stack from `caller`
	local ts
	ts="$(_polap__now)"
	local i=0
	while caller "$i" >/dev/null 2>&1; do
		local line func file
		read -r line func file < <(caller "$i")
		: "${func:=main}"
		_polap__emit 0 "[$ts ${func}@${file##*/}:${line}] #$i"
		i=$((i + 1))
	done
}

################################################################################
# ERR handler
################################################################################
_polap_err_handler() {
	local ec="${1:-1}" freeze_idx="${2:-$_POLAP_DBG_IDX}"

	# Find the most relevant command slot
	local idx
	idx="$(_polap_ring_find_last_relevant)"

	# Print the failing command line
	_polap_print_failed_cmd_at "$idx"

	# Show wrapper chain backwards from frozen head
	_polap_print_called_from_chain "$(((freeze_idx + _POLAP_DBG_MAX - 1) % _POLAP_DBG_MAX))"

	# Show a compact bash caller stack
	_polap_print_bash_caller_trace

	# Finally, emit a one-liner exit summary
	local ts
	ts="$(_polap__now)"
	# Use `caller 0` to point to the site that tripped ERR
	local line func file
	if read -r line func file < <(caller 0); then : "${func:=main}"; else
		file="${BASH_SOURCE[1]-${BASH_SOURCE[0]}}"
		line="${BASH_LINENO[0]-0}"
		func="${FUNCNAME[1]-main}"
	fi
	_polap__emit 0 "[$ts ${func}@${file##*/}:${line}] ERROR: exit ${ec}"
}

################################################################################
# Public: enable traps
################################################################################
polap_trap_enable() {
	# Inherit traps into functions/subshells
	set -o errtrace
	set -o functrace

	# Install handlers. Freeze DEBUG while printing to avoid self-noise.
	trap 'ec=$?;
        __freeze_idx="$_POLAP_DBG_IDX";
        trap - DEBUG;
        _polap_err_handler "$ec" "$__freeze_idx";
        trap "_polap_on_debug" DEBUG
      ' ERR

	trap '_polap_on_debug' DEBUG

	# Optional xtrace toggle here (leave off by default; user can enable with _POLAP_DEBUG=1)
	if [[ "${_POLAP_DEBUG}" -eq 1 ]]; then
		set -x
	fi
}
