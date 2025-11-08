#!/usr/bin/env bash
# main.sh
# Version: v0.1.0

set -Eeuo pipefail
set -o errtrace
set -o functrace

POLAPLIB="$(cd "$(dirname "${BASH_SOURCE[0]}")/polaplib" && pwd)"
source "${POLAPLIB}/polap-lib-crash.sh"
crash_enable

source "${POLAPLIB}/cmd-test.sh"

printf '[%s] Starting main test run\n' "$(_here)" >&2
# cmd_test_main r 7 "provoking nested failures across runtimes" infile.txt outfile.txt
# cmd_test_py
cmd_test_r
