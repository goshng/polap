#!/usr/bin/env bash
# polap-lib-logcallstack.sh
# Version: v0.4.0
# Usage: source "${_POLAPLIB_DIR}/polap-lib-logcallstack.sh"; polap_trap_enable

# --- Optional xtrace -----------------------------------------------------------
if [[ "${_POLAP_DEBUG:-0}" -eq 1 ]]; then
	export PS4='[${EPOCHREALTIME} ${FUNCNAME[0]}@${BASH_SOURCE##*/}:${LINENO}] '
	set -x
fi

# Inherit ERR/DEBUG into functions & subshells
set -o errtrace
set -o functrace

# --- Header + time -------------------------------------------------------------
_polap_now() { date "+%Y-%m-%d %H:%M:%S"; }

_polap_hdr() {
	# frame mapping: we’re in handler (0), caller of handler (1), failing frame’s file:line = BASH_SOURCE[1]:BASH_LINENO[0]
	local func="${FUNCNAME[2]:-main}"
	local file="$(basename "${BASH_SOURCE[2]:-${BASH_SOURCE[0]}}")"
	local line="${BASH_LINENO[1]:-0}"
	printf "[%s %s@%s:%s]" "$(_polap_now)" "$func" "$file" "$line"
}

polap_log_info() { printf "%s %s\n" "$(_polap_hdr)" "$*" >&1; }
polap_log_warn() { printf "%s WARNING: %s\n" "$(_polap_hdr)" "$*" >&2; }
polap_log_err() {
	printf "%s ERROR: %s\n" "$(_polap_hdr)" "$*" >&2
	polap_stack_dump
}

polap_log_function_start() { polap_log_info "Function start: $1"; }
polap_log_function_end() { polap_log_info "Function end:   $1"; }

# --- DEBUG ring buffer to capture the last N commands and their file:line -----
: "${_POLAP_DBG_MAX:=16}" # keep last N commands
declare -a _POLAP_DBG_CMD=() _POLAP_DBG_FILE=() _POLAP_DBG_LINE=() _POLAP_DBG_FUNC=()
_POLAP_DBG_IDX=0

_polap_on_debug() {
	# DEBUG runs *before* each simple command executes.
	# We capture the command, its source and line number.
	local idx=$_POLAP_DBG_IDX
	_POLAP_DBG_CMD[$idx]="$BASH_COMMAND"
	_POLAP_DBG_FILE[$idx]="${BASH_SOURCE[0]##*/}"
	_POLAP_DBG_LINE[$idx]="$LINENO"
	_POLAP_DBG_FUNC[$idx]="${FUNCNAME[1]:-main}"
	_POLAP_DBG_IDX=$(((idx + 1) % _POLAP_DBG_MAX))
}

# --- Stack printer: skip the duplicate caller line (#0) -----------------------
polap_stack_dump() {
	local n="${#FUNCNAME[@]}"
	# frames: 0=this func, 1=handler, 2=caller; skip the duplicate header by starting at i=3
	for ((i = 3; i < n; i++)); do
		local func="${FUNCNAME[$i]:-main}"
		local file="$(basename "${BASH_SOURCE[$i]}")"
		local line="${BASH_LINENO[$((i - 1))]:-0}"
		local depth="$((i - 3))"
		printf "[%s %s@%s:%s] #%d\n" "$(_polap_now)" "$func" "$file" "$line" "$depth" >&2
	done
}

# --- Pretty-printer for the failing command (from DEBUG ring) -----------------
_polap_print_failed_cmd() {
	# last recorded command is at index (idx - 1)
	local idx=$(((_POLAP_DBG_IDX + _POLAP_DBG_MAX - 1) % _POLAP_DBG_MAX))
	local cmd="${_POLAP_DBG_CMD[$idx]-}"
	local file="${_POLAP_DBG_FILE[$idx]-}"
	local line="${_POLAP_DBG_LINE[$idx]-}"
	local func="${_POLAP_DBG_FUNC[$idx]-}"
	if [[ -n "$cmd" && -n "$file" && -n "$line" ]]; then
		# Use %q so weird chars are visible
		printf "%s FAILED CMD at %s@%s:%s: %q\n" "$(_polap_hdr)" "${func:-main}" "${file}" "${line}" "$cmd" >&2
	else
		# Fallback: use BASH_COMMAND from ERR context
		printf "%s FAILED CMD: %q\n" "$(_polap_hdr)" "$BASH_COMMAND" >&2
	fi
}

# --- ERR trap installers -------------------------------------------------------
# v1: exit after printing (legacy behavior)
v1_polap_trap_enable() {
	trap 'ec=$?;
        _polap_print_failed_cmd;
        printf "[%s %s@%s:%s] ERROR: exit %d\n" "$(_polap_now)" "${FUNCNAME[1]:-main}" \
               "$(basename "${BASH_SOURCE[1]}")" "${BASH_LINENO[0]:-0}" "$ec" >&2;
        polap_stack_dump;
        exit "$ec"' ERR
	trap '_polap_on_debug' DEBUG
}

# Default: log + stack, do NOT exit (lets set -e decide, or callers exit)
polap_trap_enable() {
	trap 'ec=$?;
        _polap_print_failed_cmd;
        printf "[%s %s@%s:%s] ERROR: exit %d\n" "$(_polap_now)" "${FUNCNAME[1]:-main}" \
               "$(basename "${BASH_SOURCE[1]}")" "${BASH_LINENO[0]:-0}" "$ec" >&2;
        polap_stack_dump' ERR
	trap '_polap_on_debug' DEBUG
}
