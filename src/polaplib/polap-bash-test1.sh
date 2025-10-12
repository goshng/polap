#!/usr/bin/env bash
# Version: v0.1.0
# Test driver for call-stack logging and ERR trap.

set -euo pipefail
set -o errtrace
IFS=$'\n\t'

# Resolve base dir (polaplib/)
_PPOLAP="${_POLAPLIB_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)}"
export _POLAPLIB_DIR="$_PPOLAP"
export POLAPLIB_DIR="$_PPOLAP" # also export POLAPLIB_DIR for Python/R convenience

# Source the logger and enable trap
source "${_POLAPLIB_DIR}/polap-lib-logcallstack.sh"
polap_trap_enable

# Source the library under test
source "${_POLAPLIB_DIR}/polap-lib-test1.sh"

_usage() {
	cat <<EOF
Usage: $(basename "$0") [bash|python|r]
Default mode: bash

Examples:
  $(basename "$0") bash
  $(basename "$0") python
  $(basename "$0") r
EOF
}

mode="${1:-bash}"

polap_log_info "Test mode: ${mode}"

# Entry function that calls into the library and triggers the error
run_test1() {
	polap_log_info "Calling library entrypoint with mode='${mode}'"
	polap_lib_test1_entry "${mode}"
	polap_log_info "Library returned (unexpected if error path)"
}

case "$mode" in
bash | python | r) run_test1 ;;
*)
	_usage
	exit 1
	;;
esac

polap_log_info "End of test driver (no error raised)"
