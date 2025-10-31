#!/usr/bin/env bash
# Version: v0.3.1
# Usage: source "${_POLAPLIB_DIR}/polap-lib-logcallstack.sh"

# xtrace if requested
if [[ "${_POLAP_DEBUG:-0}" -eq 1 ]]; then
	export PS4='[${EPOCHREALTIME} ${FUNCNAME[0]}@${BASH_SOURCE##*/}:${LINENO}] '
	set -x
fi
set -o errtrace

_polap_now() { date "+%Y-%m-%d %H:%M:%S"; }

_polap_hdr() {
	local func="${FUNCNAME[2]:-main}"
	local file="$(basename "${BASH_SOURCE[2]:-${BASH_SOURCE[0]}}")"
	local line="${BASH_LINENO[1]:-0}"
	printf "[%s %s@%s:%s]" "$(_polap_now)" "$func" "$file" "$line"
}

polap_stack_dump() {
	local n="${#FUNCNAME[@]}"
	# frames: 0=this func, 1=polap_log_err or trap, start at 2 (caller)
	for ((i = 2; i < n; i++)); do
		local func="${FUNCNAME[$i]:-main}"
		local file="$(basename "${BASH_SOURCE[$i]}")"
		local line="${BASH_LINENO[$((i - 1))]:-0}"
		local depth="$((i - 2))"
		printf "[%s %s@%s:%s] #%d\n" "$(_polap_now)" "$func" "$file" "$line" "$depth" >&2
	done
}

polap_log_info() { printf "%s %s\n" "$(_polap_hdr)" "$*" >&1; }
polap_log_warn() { printf "%s WARNING: %s\n" "$(_polap_hdr)" "$*" >&2; }
polap_log_err() {
	printf "%s ERROR: %s\n" "$(_polap_hdr)" "$*" >&2
	polap_stack_dump
}

# optional: auto stack on uncaught ERR
v1_polap_trap_enable() {
	trap 'ec=$?;
        printf "[%s %s@%s:%s] ERROR: exit %d\n" "$(_polap_now)" "${FUNCNAME[1]:-main}" \
               "$(basename "${BASH_SOURCE[1]}")" "${BASH_LINENO[0]:-0}" "$ec" >&2;
        polap_stack_dump; exit "$ec"' ERR
}

polap_trap_enable() {
	trap 'ec=$?;
        printf "[%s %s@%s:%s] ERROR: exit %d\n" "$(_polap_now)" "${FUNCNAME[1]:-main}" \
               "$(basename "${BASH_SOURCE[1]}")" "${BASH_LINENO[0]:-0}" "$ec" >&2;
        polap_stack_dump' ERR
}
polap_log_function_start() { polap_log_info "Function start: $1"; }
polap_log_function_end() { polap_log_info "Function end:   $1"; }
