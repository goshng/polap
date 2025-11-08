#!/usr/bin/env bash
set -Eeuo pipefail
set -o errtrace
set -o functrace

# enable verbose tracing in the failsafe (prints when DEBUG/ERR fire)
export POLAP_FAILSAFE_TRACE=1

# source your current failsafe
source "${BASH_SOURCE[0]%/*}/polap-lib-failsafe.sh"
polap_enable_failsafe
echo 1

echo "[mini] next line will fail"
bash -c "exit 7"
