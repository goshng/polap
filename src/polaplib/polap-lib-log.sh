#!/usr/bin/env bash
# File: polaplib/polap-lib-log.sh
# Version: v0.10.2
################################################################################
# Contract (simple & predictable):
# • Do NOT rewire stdout/stderr globally. External tools (Bash/R/Python) print as usual.
# • Every _polap_log* API appends a timestamped line to ${LOG_FILE}.
# • Screen output obeys _arg_verbose:
#     0 → print none
#     1 → print L0 only
#     2 → print L0,L1
#     3 → print L0,L1,L2
#     4 → print L0,L1,L2,L3
# • *_cmd / *_pipe / *_pipe_no_exit only *announce* then execute as-is
#   (no global redirects). *_cmdout runs argv with stdout+stderr → ${_polap_output_dest}.
# • Enhancements over v0.10.1:
#     - Each log line includes conda env (CONDA_DEFAULT_ENV) after timestamp.
#     - *_cmd / *_pipe / *_pipe_no_exit stash POLAP_LAST_CMD_PREVIEW (+site) for ERR handler.
################################################################################

# ---------- Defaults (nounset-safe) -------------------------------------------
: "${POLAP_LOG_FILE:=./polap.log}"
: "${POLAP_VERBOSE:=1}"
: "${_POLAP_RELEASE:=0}" # 0: include [func@file:line] tag; 1: omit
# Accept legacy misspelling if present
: "${_org_verbose:=}"
if [[ -n "${_org_verbose}" && -z "${_arg_verbose:-}" ]]; then _arg_verbose="${_org_verbose}"; fi
: "${_arg_verbose:=1}"               # 0..4 screen gating
: "${POLAP_LOG_HEAD_N:=50}"          # lines for *_head
: "${_polap_output_dest:=/dev/null}" # sink for *_cmdout

# Prefer namespaced env, fall back to a sensible default
# We use LOG_FILE and POLAP_VERBOSE to control
LOG_FILE="$POLAP_LOG_FILE"

# include guard
[[ -n "${_POLAP_LOG_SOURCED:-}" ]] || _POLAP_LOG_SOURCED=1

# Ensure log file directory exists
_polap__ensure_log_dir() {
	local d
	d="$(dirname -- "$LOG_FILE")"
	[[ -d "$d" ]] || mkdir -p -- "$d" 2>/dev/null || true
}
_polap__ensure_log_dir

# ---------- Tag builder (show *caller* of public API) -------------------------
# Depth math: user → public wrapper → _emit → _tag  ⇒ depth=3
_polap__tag() {
	if [[ "${_POLAP_RELEASE}" == "1" ]]; then
		printf ""
	else
		local depth="${1:-3}"
		local func="${FUNCNAME[$depth]-main}"
		local file="$(basename -- "${BASH_SOURCE[$depth]-${BASH_SOURCE[0]}}")"
		local line="${BASH_LINENO[$((depth - 1))]-0}"
		printf "%s@%s:%s] " "$func" "$file" "$line"
	fi
}

# ---------- Screen gating by _arg_verbose -------------------------------------
# Map "Lx" → numeric 0..3
_polap__lvlnum() {
	case "$1" in
	L0) echo 0 ;;
	L1) echo 1 ;;
	L2) echo 2 ;;
	L3) echo 3 ;;
	*) echo 0 ;;
	esac
}
_polap__screen_allowed() {
	# true when _arg_verbose >= (level+1)
	local lvlnum
	lvlnum="$(_polap__lvlnum "$1")"
	local need=$((lvlnum + 1))
	[[ "${POLAP_VERBOSE:-1}" -ge "$need" ]]
}

# ---------- Timestamp + conda env ---------------------------------------------
_polap__ts_env() {
	# Example: [2025-11-06 17:18:03 polap-dev]
	local ts env
	ts="$(date '+%Y-%m-%d %H:%M:%S')"
	env="${CONDA_DEFAULT_ENV:-none}"
	printf "[%s (%s)" "$ts" "$env"
}

# ---------- Core emitter: ALWAYS log file; maybe screen per verbosity ----------
_polap__emit() {
	local lvl="${1:-L0}" depth="${2:-3}"
	shift 2 || true
	local msg="$*"
	local tag tsenv
	tag="$(_polap__tag "$depth")"
	tsenv="$(_polap__ts_env)"

	# log file (always)
	{ printf "%s %s%s\n" "$tsenv" "$tag" "$msg"; } >>"$LOG_FILE"

	# screen (stderr) if verbosity allows
	if _polap__screen_allowed "$lvl"; then
		printf "%s %s%s\n" "$tsenv" "$tag" "$msg" 1>&2
	fi
}

# ---------- Public log lines (L0..L3) -----------------------------------------
_polap_log0() { _polap__emit "L0" 3 "$*"; }
_polap_log1() { _polap__emit "L1" 3 "$*"; }
_polap_log2() { _polap__emit "L2" 3 "$*"; }
_polap_log3() { _polap__emit "L3" 3 "$*"; }

_polap_log0n() {
	local i="${1}"
	i=$((i + 3))
	_polap__emit "L0" 3 "$*"
}
_polap_log1n() {
	local i="${1}"
	i=$((i + 3))
	_polap__emit "L1" 3 "$*"
}
_polap_log2n() {
	local i="${1}"
	i=$((i + 3))
	_polap__emit "L2" 3 "$*"
}
_polap_log3n() {
	local i="${1}"
	i=$((i + 3))
	_polap__emit "L3" 3 "$*"
}

_polap_log0_log() { _polap_log0 "LOG: $*"; }
_polap_log1_log() { _polap_log1 "LOG: $*"; }
_polap_log2_log() { _polap_log2 "LOG: $*"; }
_polap_log3_log() { _polap_log3 "LOG: $*"; }

# ---------- Where-am-I helper (func:file:line) --------------------------------
_polap_here() {
	if caller 0 >/dev/null 2>&1; then
		local l f fl
		read -r l f fl < <(caller 0)
		: "${f:=main}"
		printf '%s:%s:%s' "$f" "${fl##*/}" "$l"
	else
		printf 'main:%s:%s' "${BASH_SOURCE[0]##*/}" "${LINENO:-0}"
	fi
}

# ---------- Command preview ---------------------------------------------------
_polap__preview() {
	local s="" a
	for a in "$@"; do printf -v s '%s %q' "$s" "$a"; done
	printf '%s' "${s# }"
}

# =========================
# CMD (argv) helpers
# =========================
_polap__cmd_core() {
	local lvl="$1"
	shift
	local preview
	preview="$(_polap__preview "$@")"
	_polap__emit "$lvl" 3 "command: ${preview}"

	# Stash for ERR handler to avoid 'FAILED CMD: "$@"'
	POLAP_LAST_CMD_PREVIEW="${preview}"
	POLAP_LAST_CMD_FILE="$(basename -- "${BASH_SOURCE[1]-${BASH_SOURCE[0]}}")"
	POLAP_LAST_CMD_LINE="${BASH_LINENO[0]-0}"
	export POLAP_LAST_CMD_PREVIEW POLAP_LAST_CMD_FILE POLAP_LAST_CMD_LINE

	"$@"
}
_polap_log0_cmd() { _polap__cmd_core "L0" "$@"; }
_polap_log1_cmd() { _polap__cmd_core "L1" "$@"; }
_polap_log2_cmd() { _polap__cmd_core "L2" "$@"; }
_polap_log3_cmd() { _polap__cmd_core "L3" "$@"; }

# ==============================================
# CMDOUT — run argv with stdout+stderr to DEST
# ==============================================
_polap__cmdout_core() {
	local lvl="$1"
	shift
	local dest="${_polap_output_dest:-/dev/null}"
	local preview
	preview="$(_polap__preview "$@")"
	_polap__emit "$lvl" 3 "command: ${preview} 1>${dest} 2>&1"

	# Stash preview for ERR context (even though we redirect)
	POLAP_LAST_CMD_PREVIEW="${preview} 1>${dest} 2>&1"
	POLAP_LAST_CMD_FILE="$(basename -- "${BASH_SOURCE[1]-${BASH_SOURCE[0]}}")"
	POLAP_LAST_CMD_LINE="${BASH_LINENO[0]-0}"
	export POLAP_LAST_CMD_PREVIEW POLAP_LAST_CMD_FILE POLAP_LAST_CMD_LINE

	"$@" >"${dest}" 2>&1
}
# public wrappers
_polap_log_cmdout() {
	local n="$1"
	shift
	_polap__cmdout_core "L${n}" "$@"
}
_polap_log0_cmdout() { _polap__cmdout_core "L0" "$@"; }
_polap_log1_cmdout() { _polap__cmdout_core "L1" "$@"; }
_polap_log2_cmdout() { _polap__cmdout_core "L2" "$@"; }
_polap_log3_cmdout() { _polap__cmdout_core "L3" "$@"; }

# ========================
# PIPE (string/argv) run
# ========================
_polap__pipe() {
	local lvl="$1"
	shift
	local mode="string" out_path="" append_path=""
	while [[ $# -gt 0 ]]; do
		case "$1" in
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
		--)
			shift
			break
			;;
		*) break ;;
		esac
	done
	local display
	if [[ "$mode" == "argv" ]]; then
		display="$(_polap__preview "$@")"
		[[ -n "$out_path" ]] && display+=" > ${out_path}"
		[[ -n "$append_path" ]] && display+=" >> ${append_path}"
	else
		display="$1"
	fi
	_polap__emit "$lvl" 3 "command: ${display}"

	# Stash for ERR handler
	POLAP_LAST_CMD_PREVIEW="${display}"
	POLAP_LAST_CMD_FILE="$(basename -- "${BASH_SOURCE[1]-${BASH_SOURCE[0]}}")"
	POLAP_LAST_CMD_LINE="${BASH_LINENO[0]-0}"
	export POLAP_LAST_CMD_PREVIEW POLAP_LAST_CMD_FILE POLAP_LAST_CMD_LINE

	if [[ "$mode" == "argv" ]]; then
		if [[ -n "$out_path" ]]; then
			"$@" >"$out_path"
			return $?
		elif [[ -n "$append_path" ]]; then
			"$@" >>"$append_path"
			return $?
		else
			"$@"
			return $?
		fi
	else
		eval "$1"
		return $?
	fi
}
_polap_log0_pipe() { _polap__pipe "L0" "$@"; }
_polap_log1_pipe() { _polap__pipe "L1" "$@"; }
_polap_log2_pipe() { _polap__pipe "L2" "$@"; }
_polap_log3_pipe() { _polap__pipe "L3" "$@"; }

# ==============================
# PIPE (no-exit) — tolerant run
# ==============================
_polap__pipe_no_exit() {
	local lvl="$1"
	shift
	local mode="string" out_path="" append_path="" ok_codes="" status_var=""
	while [[ $# -gt 0 ]]; do
		case "$1" in
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
		--ok | --ok-codes)
			ok_codes="$2"
			shift 2
			;;
		--status-var)
			status_var="$2"
			shift 2
			;;
		--)
			shift
			break
			;;
		*) break ;;
		esac
	done

	local display
	if [[ "$mode" == "argv" ]]; then
		display="$(_polap__preview "$@")"
		[[ -n "$out_path" ]] && display+=" > ${out_path}"
		[[ -n "$append_path" ]] && display+=" >> ${append_path}"
	else
		display="$1"
	fi
	_polap__emit "$lvl" 3 "command: ${display}"

	# Stash for ERR handler
	POLAP_LAST_CMD_PREVIEW="${display}"
	POLAP_LAST_CMD_FILE="$(basename -- "${BASH_SOURCE[1]-${BASH_SOURCE[0]}}")"
	POLAP_LAST_CMD_LINE="${BASH_LINENO[0]-0}"
	export POLAP_LAST_CMD_PREVIEW POLAP_LAST_CMD_FILE POLAP_LAST_CMD_LINE

	local had_e=0 rc=0
	case $- in *e*) had_e=1 ;; esac
	set +e
	if [[ "$mode" == "argv" ]]; then
		if [[ -n "$out_path" ]]; then
			"$@" >"$out_path"
			rc=$?
		elif [[ -n "$append_path" ]]; then
			"$@" >>"$append_path"
			rc=$?
		else
			"$@"
			rc=$?
		fi
	else
		eval "$1"
		rc=$?
	fi
	((had_e)) && set -e

	[[ -n "$status_var" ]] && printf -v "$status_var" '%s' "$rc"
	if [[ -n "$ok_codes" ]]; then
		local ok
		for ok in $ok_codes; do [[ "$rc" -eq "$ok" ]] && return 0; done
	fi
	return "$rc"
}
_polap_log0_pipe_no_exit() { _polap__pipe_no_exit "L0" "$@"; }
_polap_log1_pipe_no_exit() { _polap__pipe_no_exit "L1" "$@"; }
_polap_log2_pipe_no_exit() { _polap__pipe_no_exit "L2" "$@"; }
_polap_log3_pipe_no_exit() { _polap__pipe_no_exit "L3" "$@"; }

# ---------- File helpers (CAT / HEAD / COLUMN) --------------------------------
_polap__file_info() {
	local path="$1" label="$2"
	if [[ -f "$path" ]]; then
		local n
		n=$(wc -l <"$path" 2>/dev/null || echo 0)
		_polap__emit "$label" 3 "FILE (${n} lines): ${path}"
		return 0
	else
		_polap__emit "$label" 3 "ERROR: no such file: ${path}"
		return 1
	fi
}
_polap_log0_cat() { _polap__file_info "${1:?}" "L0" && cat -- "$1" 1>&2 || true; }
_polap_log1_cat() { _polap__file_info "${1:?}" "L1" && cat -- "$1" 1>&2 || true; }
_polap_log2_cat() { _polap__file_info "${1:?}" "L2" && cat -- "$1" 1>&2 || true; }
_polap_log3_cat() { _polap__file_info "${1:?}" "L3" && cat -- "$1" 1>&2 || true; }

_polap__head_n() {
	local n="${POLAP_LOG_HEAD_N:-50}"
	[[ "$n" =~ ^[0-9]+$ ]] || n=50
	printf '%s' "$n"
}
_polap_log0_head() { _polap__file_info "${1:?}" "L0" && head -n "$(_polap__head_n)" -- "$1" 1>&2 || true; }
_polap_log1_head() { _polap__file_info "${1:?}" "L1" && head -n "$(_polap__head_n)" -- "$1" 1>&2 || true; }
_polap_log2_head() { _polap__file_info "${1:?}" "L2" && head -n "$(_polap__head_n)" -- "$1" 1>&2 || true; }
_polap_log3_head() { _polap__file_info "${1:?}" "L3" && head -n "$(_polap__head_n)" -- "$1" 1>&2 || true; }

_polap_log0_column() { _polap__file_info "${1:?}" "L0" && column -t -- "$1" 1>&2 || true; }
_polap_log1_column() { _polap__file_info "${1:?}" "L1" && column -t -- "$1" 1>&2 || true; }
_polap_log2_column() { _polap__file_info "${1:?}" "L2" && column -t -- "$1" 1>&2 || true; }
_polap_log3_column() { _polap__file_info "${1:?}" "L3" && column -t -- "$1" 1>&2 || true; }

# ---------- Function banner ---------------------------------------------------
_polap_log_function() { _polap__emit "L3" 3 "$*"; }

_polap_log0_file() {
	_polap_log0 "FILE: $@"
}

_polap_log1_file() {
	_polap_log1 "FILE: $@"
}

_polap_log2_file() {
	_polap_log2 "FILE: $@"
}

_polap_log3_file() {
	_polap_log3 "FILE: $@"
}
