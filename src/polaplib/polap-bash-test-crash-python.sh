#!/usr/bin/env bash
# polap-bash-test-crash-python.sh
# Version: v0.1.0
set -Eeuo pipefail
: "${_POLAPLIB_DIR:=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
source "${_POLAPLIB_DIR}/polap-lib-logcallstack.sh"
polap_trap_enable

# Call Python directly; no PYTHONPATH/sitecustomize needed
python3 "${_POLAPLIB_DIR}/polap-py-test-crash.py" --msg "from bash wrapper" --code 99
