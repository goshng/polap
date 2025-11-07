#!/usr/bin/env bash
# polaplib/polap-lib-log.sh
# Version: v0.9.0
################################################################################
# This file is part of polap.
#
# Goal (simplified & predictable):
#   • Do NOT rewire stdout/stderr globally. Let Bash/R/Python print as usual.
#   • _polap_log_* helpers write the same line to:
#       - stderr (so you see it on screen), and
#       - ${LOG_FILE} (append, with timestamp)
#   • Command wrappers (_polap_log*_cmd) only *announce* the command, then run
#     it normally so its own stdout/stderr behave exactly as usual.
#
# Assumptions:
#   • polap-variables-main.sh sets LOG_FILE (preferred).
#     If not set, we default to "./polap.log".
################################################################################

# ----- Defaults (nounset-safe) ------------------------------------------------
: "${LOG_FILE:=./polap.log}"
: "${_POLAP_RELEASE:=0}"   # 0 = dev (tag func@file:line), 1 = release (no tag)
: "${_arg_verbose:=1}"     # 0..4 (only affects *screen* gate; log file always appends)
: "${_arg_log_stderr:=on}" # ignored here; we always write to stderr + file

# ----- Internal: build a tag (dev mode prints func@file:line) -----------------
_polap__tag() {
	if [[ "${_POLAP_RELEASE}" == "1" ]]; then
		printf ""
	else
		# caller is the function that invoked the logger (depth=2)
		local depth="${1:-2}"
		local func="${FUNCNAME[$depth]-main}"
		local file="$(basename -- "${BASH_SOURCE[$depth]-${BASH_SOURCE[0]}}")"
		local line="${BASH_LINENO[$((depth - 1))]-0}"
		printf "[%s@%s:%s] " "$func" "$file" "$line"
	fi
}

# ----- Internal: write one line to stderr and to LOG_FILE ---------------------
# Args:
#   $1 = level label (e.g., "L0" or "L3") – informational only
#   $2 = tag depth (defaults to 2)
#   $3... = message
_polap__emit() {
	local lvl="${1:-L0}" depth="${2:-2}"
	shift 2 || true
	local msg="$*"
	local tag
	tag="$(_polap__tag "$depth")"
	local ts
	ts="$(date '+%Y-%m-%d %H:%M:%S')"

	# stderr (screen)
	printf "%s%s\n" "$tag" "$msg" 1>&2

	# log file (append with timestamp)
	# shellcheck disable=SC2129
	{
		printf "[%s] %s%s\n" "$ts" "$tag" "$msg"
	} >>"$LOG_FILE"
}

# ----- Public log lines (visible to screen and file) --------------------------
_polap_log0() { _polap__emit "L0" 2 "$*"; } # minimal / user-facing
_polap_log1() { _polap__emit "L1" 2 "$*"; } # info
_polap_log2() { _polap__emit "L2" 2 "$*"; } # detail
_polap_log3() { _polap__emit "L3" 2 "$*"; } # debug

# Back-compat short aliases used elsewhere (write both places)
_polap_log0_log() { _polap_log0 "LOG: $*"; }
_polap_log1_log() { _polap_log1 "LOG: $*"; }
_polap_log2_log() { _polap_log2 "LOG: $*"; }
_polap_log3_log() { _polap_log3 "LOG: $*"; }

# ----- Compact "where am I" helper for your headers (func:file:line) ----------
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

# ----- Command preview (no redirection; run as-is) ----------------------------
# Prints a quoted, single-line preview to stderr and to LOG_FILE, then executes.
# Usage: _polap_logN_cmd cmd arg1 arg2 ...
_polap__preview() {
	local s="" a
	for a in "$@"; do printf -v s '%s %q' "$s" "$a"; done
	printf "%s" "${s# }"
}

_polap__cmd_core() {
	local lvl="$1"
	shift
	local preview
	preview="$(_polap__preview "$@")"
	_polap__emit "$lvl" 3 "command: ${preview}"
	"$@"
}

_polap_log0_cmd() { _polap__cmd_core "L0" "$@"; }
_polap_log1_cmd() { _polap__cmd_core "L1" "$@"; }
_polap_log2_cmd() { _polap__cmd_core "L2" "$@"; }
_polap_log3_cmd() { _polap__cmd_core "L3" "$@"; }

# ----- No-exit variant (returns real status but does not kill caller) ---------
# Usage:
#   _polap_log*_pipe_no_exit --argv cmd arg1 ...         # argv-safe
#   _polap_log*_pipe_no_exit "string to eval"            # eval mode
_polap__pipe_no_exit() {
	local lvl="$1"
	shift
	local mode="string" status_var="" out_path="" append_path="" ok_codes=""
	while [[ $# -gt 0 ]]; do
		case "$1" in
		--argv)
			mode="argv"
			shift
			;;
		--status-var)
			status_var="$2"
			shift 2
			;;
		--ok | --ok-codes)
			ok_codes="$2"
			shift 2
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

	local st=0 had_e=0
	case $- in *e*) had_e=1 ;; esac
	set +e

	if [[ "$mode" == "argv" ]]; then
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

	((had_e)) && set -e
	[[ -n "$status_var" ]] && printf -v "$status_var" '%s' "$st"

	if [[ -n "$ok_codes" ]]; then
		for ok in $ok_codes; do [[ "$st" -eq "$ok" ]] && return 0; done
	fi
	return "$st"
}

_polap_log0_pipe_no_exit() { _polap__pipe_no_exit "L0" "$@"; }
_polap_log1_pipe_no_exit() { _polap__pipe_no_exit "L1" "$@"; }
_polap_log2_pipe_no_exit() { _polap__pipe_no_exit "L2" "$@"; }
_polap_log3_pipe_no_exit() { _polap__pipe_no_exit "L3" "$@"; }

# ----- File helpers (screen + LOG_FILE; never redirect file contents globally)-
_polap_log*_common_cat() {
	local lvl="$1" path="$2" hdr="$3"
	if [[ -f "$path" ]]; then
		local n
		n=$(wc -l <"$path" 2>/dev/null || echo 0)
		_polap__emit "$lvl" 2 "${hdr} (${n} lines): $path"
		# Show head few lines to screen; always append full notice to log already
		# (Do not cat to file; we do not mirror full file contents into LOG_FILE)
		head -n 50 -- "$path" 1>&2 || true
	else
		_polap__emit "$lvl" 2 "ERROR: no such file: $path"
	fi
}

_polap_log0_cat() { _polap_log*_common_cat "L0" "$1" "FILE"; }
_polap_log1_cat() { _polap_log*_common_cat "L1" "$1" "FILE"; }
_polap_log2_cat() { _polap_log*_common_cat "L2" "$1" "FILE"; }
_polap_log3_cat() { _polap_log*_common_cat "L3" "$1" "FILE"; }

_polap_log0_head() {
	local lvl="L0" path="$1"
	if [[ -f "$path" ]]; then
		local n
		n=$(wc -l <"$path" 2>/dev/null || echo 0)
		_polap__emit "$lvl" 2 "FILE head (${n} lines): $path"
		head -- "$path" 1>&2 || true
	else
		_polap__emit "$lvl" 2 "ERROR: no such file: $path"
	fi
}
_polap_log1_head() { _polap_log0_head "$@"; }
_polap_log2_head() { _polap_log0_head "$@"; }
_polap_log3_head() { _polap_log0_head "$@"; }

_polap_log0_file() { _polap_log0 "FILE: $*"; }
_polap_log1_file() { _polap_log1 "FILE: $*"; }
_polap_log2_file() { _polap_log2 "FILE: $*"; }
_polap_log3_file() { _polap_log3 "FILE: $*"; }

# ----- Function banner (start/end markers) -----------------------------------
_polap_log_function() { _polap__emit "L1" 2 "$*"; }

################################################################################
# Notes:
# • No FD juggling here. stdout/stderr for all tools (R/Python/Bash) keep their
#   natural behavior. Only explicit _polap_log_* calls are duplicated to stderr
#   and ${LOG_FILE}.
# • If you want timestamps on every external tool line, prefer enabling those
#   tools’ own verbose/timestamp flags rather than global redirection here.
################################################################################
