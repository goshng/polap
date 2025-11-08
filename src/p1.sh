#!/usr/bin/env bash

set -Eeuo pipefail
set -o errtrace
set -o functrace

# Must be executed (not sourced)
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
	echo "[ERROR] This script must be executed, not sourced: use 'bash $BASH_SOURCE'" >&2
	return 1 2>/dev/null || exit 1
fi

: "${_POLAP_DEBUG:=0}"
export _POLAP_DEBUG
: "${_POLAP_RELEASE:=0}"
export _POLAP_RELEASE

POLAPLIB="$(cd "$(dirname "${BASH_SOURCE[0]}")/polaplib" && pwd)"
_POLAPLIB_DIR="${POLAPLIB}"

# source "${POLAPLIB}/polap-lib-crash.sh"
source "${_POLAPLIB_DIR}/polap-lib-version.sh"
source "${_POLAPLIB_DIR}/polap-parsing.sh"

# Autoload libraries + traps
# Do not shift "$@" here—just detect the flag.
for __a in "$@"; do
	[[ "$__a" == "--trace-load" ]] && export POLAP_AUTOLOAD_TRACE=1
done
unset __a

source "${_POLAPLIB_DIR}/polap-bash-autoload.sh"

polap_autoload "${_POLAPLIB_DIR}"

# polap_trap_enable # installs DEBUG/ERR ring printer
# crash_enable
polap_enable_failsafe

# ───────────────────────────────────────────────────────────────────────────────
# MAIN
# ───────────────────────────────────────────────────────────────────────────────

# Show help if no arguments
# if [[ $# -eq 0 ]]; then
# 	print_help
# 	touch make-menus
# 	exit "${EXIT_SUCCESS:-0}"
# fi

if [[ -z "${BASH_VERSION:-}" ]]; then
	echo "ERROR: this script must be run with bash" >&2
	exit 2
fi

# Prepare output locations and the log file path
source "${_POLAPLIB_DIR}/polap-variables-main.sh"
_polap_timer_reset
_polap_lib_log-init # sets $LOG_FILE and ensures ${_arg_outdir}/{tmp,log} exist

CMD="$0 $*"
{
	echo "POLAP: ${_polap_version}"
	echo "CMD: $CMD"
} >>"$LOG_FILE"

# Print all parsed globals from polap-parsing.sh to screen and log
# (kept minimal; still useful when debugging parsing)
# {
# 	set +u
# 	for var in $(compgen -v _arg_); do
# 		echo "$var=${!var}"
# 	done
# 	set -u
# } >>"$LOG_FILE"

# Postprocess command-line arguments
# Normalize options that need preprocessing
if [[ -n "${_arg_genomesize:-}" ]]; then
	_arg_genomesize=$(_polap_lib_unit-convert_to_int "${_arg_genomesize}")
fi

_polap_debug_log0 "verbose level: ${_arg_verbose}"

# If you want log functions to prefer FD3 instead of STDERR, ensure this is "off"
# so release-mode _polap_log* send to FD3 (the log file) rather than FD2.
: "${_arg_log_stderr:=off}"

# Call polap subcommand
"_run_polap_${_arg_menu[0]}"

[[ "${_arg_clock:-off}" == "on" ]] && date +"%Y-%m-%d %H:%M:%S" | tee -a "$LOG_FILE" >&1

echo "$(_polap_timer_log)" >>"$LOG_FILE"
