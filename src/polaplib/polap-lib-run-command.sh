#!/usr/bin/env bash
# File: polap-lib-run-command.sh
# Version: v0.2.5
# Purpose:
#   Unified Polap step runner (idempotent, dry-run, logged) with safe backend execution.
#   Delegates logging to polap-lib-log.sh (_polap_log0..3).
#   Auto-sources argv/simple backends when available.

# ── include guard ─────────────────────────────────────────────────────────────
if [[ -n "${_POLAP_RUN_COMMAND_SOURCED:-}" ]]; then
	return 0
fi
_Polap_Run_Command_Version="v0.2.5"
readonly _Polap_Run_Command_Version
_Polap_Run_Command_File="${BASH_SOURCE[0]:-polap-lib-run-command.sh}"
readonly _Polap_Run_Command_File
_POLAP_RUN_COMMAND_SOURCED=1

# ── locate Polap lib dir ──────────────────────────────────────────────────────
_POLAPLIB_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"

# ── source core logging ───────────────────────────────────────────────────────
# The log lib defines: _polap_log0 .. _polap_log3 and their control knobs.
if [[ -r "${_POLAPLIB_DIR}/polap-lib-log.sh" ]]; then
	# shellcheck source=polap-lib-log.sh
	source "${_POLAPLIB_DIR}/polap-lib-log.sh"
else
	echo "[FATAL] Missing required log library: ${_POLAPLIB_DIR}/polap-lib-log.sh" >&2
	exit 2
fi

# ── config knobs (runner-specific) ────────────────────────────────────────────
: "${POLAP_FORCE:=0}"          # 1: ignore step OK markers and re-run
: "${POLAP_DRYRUN:=0}"         # 1: print what would run, don't run
: "${POLAP_RUN_BACKEND:=auto}" # auto|argv|simple|direct

# ── tiny helpers ──────────────────────────────────────────────────────────────
polap__ensure_dir() { install -d -- "$1"; }
polap__mkparent() { install -d -- "$(dirname -- "$1")"; }

polap_die() {
	_polap_log0 "ERROR: $*"
	exit 2
}
polap_require() {
	local x
	for x in "$@"; do command -v "$x" >/dev/null 2>&1 || polap_die "Missing dependency: $x"; done
}

# ── lock helpers ──────────────────────────────────────────────────────────────
polap_lock_acquire_wrapper() {
	local lockdir="$1"
	[[ -z "$lockdir" ]] && polap_die "lock path required"
	if mkdir "$lockdir" 2>/dev/null; then
		printf 'host=%s\npid=%s\nts=%s\n' "$(hostname)" "$$" "$(date '+%Y-%m-%d %H:%M:%S')" >"$lockdir/owner"
		trap 'polap_lock_release_wrapper "'"$lockdir"'"' EXIT INT TERM
		_polap_log1 "LOCK acquired: $lockdir"
		return 0
	fi
	if [[ -s "$lockdir/owner" ]]; then
		_polap_log0 "LOCK busy ($lockdir): $(tr '\n' ' ' <"$lockdir/owner")"
	else
		_polap_log0 "LOCK busy: $lockdir"
	fi
	return 3
}

polap_lock_release_wrapper() {
	local lockdir="$1"
	if [[ -z "$lockdir" ]]; then
		_polap_log2 "LOCK release: empty path (nothing to do)"
		return 0
	fi
	if [[ -d "$lockdir" ]]; then
		if rm -rf -- "$lockdir"; then
			_polap_log1 "LOCK released: $lockdir"
		else
			_polap_log0 "WARN: failed to remove lockdir: $lockdir"
		fi
	else
		_polap_log2 "LOCK not present: $lockdir"
	fi
	return 0
}

# ── source run backends if available ──────────────────────────────────────────
# [[ -r "${_POLAPLIB_DIR}/polap-lib-run-simple.sh" ]] &&
source "${_POLAPLIB_DIR}/polap-lib-run-simple.sh"
# || true
# [[ -r "${_POLAPLIB_DIR}/polap-lib-run-argv.sh" ]] &&
source "${_POLAPLIB_DIR}/polap-lib-run-argv.sh"
# || true

# # ── safe fallback backend definitions ─────────────────────────────────────────
# if ! declare -F polap_run_argv >/dev/null 2>&1; then
# 	polap_run_argv() {
# 		local tag=""
# 		while [[ $# -gt 0 ]]; do
# 			case "$1" in
# 			--tag)
# 				tag="$2"
# 				shift 2
# 				;;
# 			--)
# 				shift
# 				break
# 				;;
# 			*) break ;;
# 			esac
# 		done
# 		local preview
# 		preview="$(printf '%q ' "$@")"
# 		preview="${preview# }"
# 		_polap_log1 "[RUN:${tag}] ${preview}"
# 		"$@"
# 	}
# fi
#
# if ! declare -F polap_run_simple >/dev/null 2>&1; then
# 	polap_run_simple() {
# 		local tag="" cmd=""
# 		while [[ $# -gt 0 ]]; do
# 			case "$1" in
# 			--tag)
# 				tag="$2"
# 				shift 2
# 				;;
# 			--)
# 				shift
# 				break
# 				;;
# 			*) break ;;
# 			esac
# 		done
# 		cmd="$*"
# 		_polap_log1 "[RUN:${tag}] ${cmd}"
# 		eval "$cmd"
# 	}
# fi

# ── backend selector ──────────────────────────────────────────────────────────
# ── backend selector ──────────────────────────────────────────────────────────
v1_polap__exec_with_backend() {
	local backend="${POLAP_RUN_BACKEND:-auto}"
	local tag="$1"
	shift
	case "$backend" in
	argv)
		polap_run_argv --tag "$tag" -- "$@"
		;;
	simple)
		# IMPORTANT: pass the *raw* string (no %q) so eval sees real pipes/redirs
		local cmd="$*"
		polap_run_simple --tag "$tag" -- "$cmd"
		;;
	direct)
		"$@"
		;;
	auto | *)
		local joined="$*"
		# Heuristically detect shell metacharacters that require eval parsing
		if [[ "$joined" =~ [\|\>\<\&\;\(\)] ]]; then
			# contains | > < & ; ( ) → run via simple backend
			if declare -F polap_run_simple >/dev/null 2>&1; then
				_polap_log2 "auto-detected shell syntax in '$tag': using simple backend"
				polap_run_simple --tag "$tag" -- "$joined"
			else
				_polap_log2 "shell syntax detected but no simple backend; falling back to bash -c"
				bash -c "$joined"
			fi
		else
			# no shell syntax → use argv backend (safest)
			if declare -F polap_run_argv >/dev/null 2>&1; then
				polap_run_argv --tag "$tag" -- "$@"
			else
				"$@"
			fi
		fi
		;;
	esac
}

v1_polap__exec_with_backend() {
	local backend="${POLAP_RUN_BACKEND:-auto}"
	local tag="$1"
	shift
	case "$backend" in
	argv)
		polap_run_argv --tag "$tag" -- "$@"
		;;
	simple)
		local cmd="$*"
		polap_run_simple --tag "$tag" -- "$cmd"
		;;
	direct)
		"$@"
		;;
	auto | *)
		if (($# == 1)); then
			# One string → the user intends shell parsing (pipes/redirects/etc.)
			polap_run_simple --tag "$tag" -- "$1"
		else
			# Multiple argv tokens → safest path, literal '|' stays data
			if declare -F polap_run_argv >/dev/null 2>&1; then
				polap_run_argv --tag "$tag" -- "$@"
			else
				"$@"
			fi
		fi
		;;
	esac
}

# ── Step-aware wrapper (idempotent, logged, dry-run) ──────────────────────────
v1_polap_run_wrapper() {
	local name="$1"
	shift
	local stepdir="${OUTDIR:?OUTDIR not set}/.polap-run/$name"
	local steplog="$stepdir/step.log"
	polap__ensure_dir "$stepdir"

	# Skip completed step unless forced
	if ((POLAP_FORCE == 0)) && [[ -e "$stepdir/.ok" ]]; then
		_polap_log1 "SKIP (ok): $name"
		return 0
	fi

	printf '%q ' "$@" >"$stepdir/cmd.txt"
	: >"$stepdir/.running"
	_polap_log1 "RUN : $name"

	if ((POLAP_DRYRUN == 1)); then
		rm -f "$stepdir/.running"
		return 0
	fi

	local rc=0
	(
		set -Eeuo pipefail
		polap__exec_with_backend "$name" "$@"
	) >>"$steplog" 2>&1 || rc=$?

	if ((rc != 0)); then
		rm -f "$stepdir/.running"
		_polap_log0 "FAIL: $name (rc=$rc) see $steplog"
		return "$rc"
	fi

	touch "$stepdir/.ok"
	rm -f "$stepdir/.running"
	_polap_log1 "DONE: $name"
}

# ---------------------------------------------------------------------------
# Single entrypoint: step orchestration + AUTO execution selection
# Assumes: polap_run_simple() and polap_run_argv() are AVAILABLE.
# Usage:
#   polap_run_wrapper <step-name> <argv...>          # argv-safe (no shell parsing)
#   polap_run_wrapper <step-name> "<shell string>"   # shell string (pipes/redirs)
# Env:
#   OUTDIR (required), POLAP_FORCE=0/1, POLAP_DRYRUN=0/1
# Creates:
#   $OUTDIR/.polap-run/<step>/{cmd.txt,step.log,.running,.ok}
# Behavior:
#   - Idempotent via .ok (skip unless POLAP_FORCE=1)
#   - Dry-run via POLAP_DRYRUN=1 (no payload execution)
#   - Strict execution in subshell (set -Eeuo pipefail)
#   - AUTO rule:
#       * if exactly 1 token after <step-name>  → run as shell string (polap_run_simple)
#       * else                                  → run argv-safe (polap_run_argv)
# ---------------------------------------------------------------------------
polap_run_wrapper() {
	local name="$1"
	shift
	local stepdir="${OUTDIR:?OUTDIR not set}/.polap-run/$name"
	local steplog="$stepdir/step.log"
	polap__ensure_dir "$stepdir"

	# Idempotency gate
	if ((POLAP_FORCE == 0)) && [[ -e "$stepdir/.ok" ]]; then
		_polap_log1n 1 "SKIP (ok): $name"
		return 0
	fi

	# Record exact argv for reproducibility
	printf '%q ' "$@" >"$stepdir/cmd.txt"
	: >"$stepdir/.running"
	_polap_log1n 1 "RUN : $name"

	# Dry-run: announce only
	if ((POLAP_DRYRUN == 1)); then
		rm -f "$stepdir/.running"
		return 0
	fi

	local rc=0
	(
		set -Eeuo pipefail

		if (($# == 1)); then
			# One token → caller intends a shell command string (pipes/redirs/etc.)
			polap_run_simple --tag "$name" -- "$1"
		else
			# Multiple tokens → argv-safe; '|' inside args remains literal data
			polap_run_argv --tag "$name" -- "$@"
		fi
	) >>"$steplog" 2>&1 || rc=$?

	if ((rc != 0)); then
		rm -f "$stepdir/.running"
		_polap_log0n 1 "FAIL: $name (rc=$rc) see $steplog"
		return "$rc"
	fi

	touch "$stepdir/.ok"
	rm -f "$stepdir/.running"
	_polap_log1n 1 "DONE: $name"
	return 0
}
