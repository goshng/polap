#!/usr/bin/env bash
# polap-bash-oatk-log.sh
# Minimal logging helpers for POLAP scripts.
#
# Features:
#   • Verbosity levels: 0=quiet, 1=info, 2=trace
#   • Logs include file:function:line
#   • Command runner logs full, safely quoted command + exit code
#   • Optional xtrace at level >=2 with readable PS4
#
# Usage (in your script):
#   source "${_POLAPLIB_DIR}/polap-bash-oatk-log.sh"
#   OATK_LOG_LEVEL=1      # or call: oatk_log_set_level 1
#   oatk_log 1 "hello"
#   oatk_enable_trace_if_needed
#   oatk_run some_cmd --flag "value"
#
# Optional: parse -v/--verbose occurrences to bump level
#   vcount=$(printf '%s\n' "$@" | grep -Eo '(^|[[:space:]])(-v|--verbose)([[:space:]]|$)' | wc -l)
#   oatk_log_set_level "$vcount"

# ---------------- defaults ----------------
: "${OATK_LOG_LEVEL:=0}" # 0=quiet, 1=info, 2=trace

# ---------------- API ----------------

oatk_log_set_level() {
	# $1: integer level (0/1/2)
	local lvl="${1:-0}"
	if [[ "$lvl" =~ ^[0-9]+$ ]]; then
		OATK_LOG_LEVEL="$lvl"
	fi
}

oatk_log_inc_level() {
	OATK_LOG_LEVEL=$((OATK_LOG_LEVEL + 1))
}

oatk_enable_trace_if_needed() {
	# Enable xtrace with clear PS4 when level >= 2
	if ((OATK_LOG_LEVEL >= 2)); then
		# file:function:line -> command
		export PS4='+${BASH_SOURCE##*/}:${FUNCNAME[0]:-main}:${LINENO}: '
		set -x
	fi
}

oatk_log() {
	# $1: level (int), $@ message
	local lvl="${1:-0}"
	shift || true
	((OATK_LOG_LEVEL >= lvl)) || return 0
	# caller info: [file:function:line] message
	# index 1 = direct caller of oatk_log
	local src="${BASH_SOURCE[1]##*/}"
	local fn="${FUNCNAME[1]:-main}"
	local ln="${BASH_LINENO[0]:-0}"
	printf '[%s:%s:%s] %s\n' "$src" "$fn" "$ln" "$*" >&2
}

oatk_cmdstr() {
	# echo a safely shell-quoted command line from $@
	local out=() a
	for a in "$@"; do
		printf -v a '%q' "$a"
		out+=("$a")
	done
	printf '%s' "${out[*]}"
}

# oatk_run() {
# 	# log command (level 1), run it, log exit code; return child status
# 	local cmdline
# 	cmdline="$(oatk_cmdstr "$@")"
# 	oatk_log 1 "exec: ${cmdline}"
# 	"$@"
# 	local st=$?
# 	oatk_log 1 "exit(${st}): ${cmdline}"
# 	return "$st"
# }
#
oatk_run() {
	# Build a safely-quoted command line
	local cmdline
	cmdline="$(oatk_cmdstr "$@")"

	# Compute caller context (file:function:line) of the site that invoked oatk_run
	# Stack inside this function:
	#   BASH_SOURCE[0] = polap-bash-oatk-log.sh (this file)
	#   BASH_SOURCE[1] = caller script file      (what we want to show)
	#   FUNCNAME[1]    = caller function name
	#   BASH_LINENO[0] = line number in caller where oatk_run was called
	local src="${BASH_SOURCE[1]##*/}"
	local fn="${FUNCNAME[1]:-main}"
	local ln="${BASH_LINENO[0]:-0}"

	# Log "exec" at info level regardless of OATK_LOG_LEVEL gate in oatk_log
	if ((OATK_LOG_LEVEL >= 1)); then
		printf '[%s:%s:%s] exec: %s\n' "$src" "$fn" "$ln" "$cmdline" >&2
	fi

	# Run the command
	"$@"
	local st=$?

	# Log exit with the same caller context
	if ((OATK_LOG_LEVEL >= 1)); then
		printf '[%s:%s:%s] exit(%d): %s\n' "$src" "$fn" "$ln" "$st" "$cmdline" >&2
	fi

	return "$st"
}

oatk_provenance_line() {
	# Print a single-line provenance summary:
	#   config=<path> overrides=<args...>
	# $1: resolved config path (or '?')
	# $@: override args to record (already filtered if needed)
	local cfg="${1:-?}"
	shift || true
	printf '[provenance] config=%s overrides=%s\n' "$cfg" "${*:-'(none)'}"
}

oatk_provenance_append() {
	# Append provenance block to a file (create dir if needed)
	# $1: file path to append
	# $2: resolved config path
	# $@: override args to record
	local fp="${1:-}"
	shift || true
	local cfg="${1:-?}"
	shift || true
	[[ -n "$fp" ]] || return 0
	mkdir -p "$(dirname -- "$fp")" 2>/dev/null || true
	{
		printf '%s  config=%s\n' "$(date -Iseconds)" "$cfg"
		printf 'overrides: %s\n\n' "${*:-'(none)'}"
	} >>"$fp"
}

# -------- optional helpers for boolean normalization --------
oatk_bool() {
	# normalize to literal "true"/"false" (leave unknowns unchanged)
	local v="${1:-}"
	case "${v,,}" in
	true | 1 | yes | y | on) echo "true" ;;
	false | 0 | no | n | off | "") echo "false" ;;
	*) echo "$v" ;;
	esac
}

# -------- optional helper to build filtered override list --------
oatk_filter_overrides() {
	# Prints "$@" but removes addressing/verbosity flags often present in POLAP launchers.
	# Use like: filtered=$(oatk_filter_overrides "${ORIG_ARGS[@]}")
	local out=() a
	for a in "$@"; do
		case "$a" in
		--path | --config-dir | --preset | -v | --verbose) continue ;;
		*) out+=("$a") ;;
		esac
	done
	printf '%s\n' "${out[@]}"
}
