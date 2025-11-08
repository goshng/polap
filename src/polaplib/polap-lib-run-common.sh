#!/usr/bin/env bash
# polaplib/polap-lib-run-common.sh
# Minimal shared utilities for polap_* run/step helpers.

# ── include guard ─────────────────────────────────────────────────────────────
if [[ -n "${_POLAP_RUN_COMMON_SOURCED:-}" ]]; then
	return 0
fi
_Polap_Run_Common_Version="v0.1.0"
readonly _Polap_Run_Common_Version
_Polap_Run_Common_File="${BASH_SOURCE[0]:-polap-lib-run-common.sh}"
readonly _Polap_Run_Common_File
_POLAP_RUN_COMMON_SOURCED=1

# ── config knobs ──────────────────────────────────────────────────────────────
: "${POLAP_VERBOSE:=1}" # 0..4 → which levels go to screen (0=quiet)
: "${LOG_FILE:=}"       # if set, all messages also append to this file
: "${POLAP_TAGDATE:=1}" # show date in logs (1=yes)

# ── small helpers ─────────────────────────────────────────────────────────────
polap__ts() { date '+%Y-%m-%d %H:%M:%S'; }
polap__envtag() {
	local e="${CONDA_DEFAULT_ENV-}"
	[[ -n "$e" ]] && printf '(%s) ' "$e"
}
# tag for the *caller* at stack depth (default = 2: your wrapper calls this)
polap__tag() {
	local depth="${1:-2}"
	local func="${FUNCNAME[$depth]-main}"
	local file="$(basename -- "${BASH_SOURCE[$depth]-${BASH_SOURCE[0]}}")"
	local line="${BASH_LINENO[$((depth - 1))]-0}"
	printf '%s@%s:%s' "$func" "$file" "$line"
}

# join/quote argv for preview (no eval)
polap__preview_argv() {
	local s="" a
	for a in "$@"; do printf -v s '%s %q' "$s" "$a"; done
	printf '%s' "${s# }"
}

# "expand" a string without executing it: only parameter expansion
# (Note: this does *not* run commands; it just expands $vars within the current shell)
polap__expand_string() {
	local s="$*"
	# shellcheck disable=SC2016
	eval "printf '%s' \"$s\""
}

# core emitter: L0..L3 → writes to LOG_FILE (if set) and maybe to screen
# depth: which stack level to print as caller tag (3 works well for wrappers)
polap__emit() {
	local level="$1" depth="$2"
	shift 2
	local msg="$*"
	local ts=""
	[[ "$POLAP_TAGDATE" == "1" ]] && ts="$(polap__ts)"
	local tag="$(polap__tag "$depth")"
	local env="$(polap__envtag)"
	local line
	if [[ -n "$ts" ]]; then
		printf -v line '[%s %s%s] %s' "$ts" "$env" "$tag" "$msg"
	else
		printf -v line '[%s%s] %s' "$env" "$tag" "$msg"
	fi

	# To file
	if [[ -n "$LOG_FILE" ]]; then
		printf '%s\n' "$line" >>"$LOG_FILE"
	fi
	# To screen depending on POLAP_VERBOSE vs level (L0=need 1, L1=2, L2=3, L3=4)
	local need
	case "$level" in
	L0) need=1 ;;
	L1) need=2 ;;
	L2) need=3 ;;
	L3) need=4 ;;
	*) need=2 ;;
	esac
	if ((POLAP_VERBOSE >= need)); then
		printf '%s\n' "$line" >&2
	fi
}

# public log shorthands (depth=3 is suited for wrappers that call emit)
polap_log0() { polap__emit L0 3 "$*"; }
polap_log1() { polap__emit L1 3 "$*"; }
polap_log2() { polap__emit L2 3 "$*"; }
polap_log3() { polap__emit L3 3 "$*"; }

# tiny file helpers
polap__ensure_dir() { install -d -- "$1"; }
polap__mkparent() { install -d -- "$(dirname -- "$1")"; }

# assert helpers (exit 1 on failure)
polap_assert_file() {
	local f="$1"
	[[ -s "$f" ]] || {
		polap__emit L0 3 "ASSERT file missing: $f"
		return 1
	}
}
polap_assert_contains() {
	local needle="$1" hay="$2"
	grep -Fq -- "$needle" "$hay" || {
		polap__emit L0 3 "ASSERT text not found: $needle in $hay"
		return 1
	}
}
