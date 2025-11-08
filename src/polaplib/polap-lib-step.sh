#!/usr/bin/env bash
# polaplib/polap-lib-step.sh
# polap_step: run a *string* command as a named step.

# shellcheck disable=SC1091
source "${BASH_SOURCE[0]%/*}/polap-lib-run-common.sh"

# Usage: polap_step "<step-name>" "<cmd string>"
polap_step() {
	local name="$1" cmd="$2"
	local start end rc
	start=$(date +%s)
	polap_log1 "▶ step: ${name}"
	polap_log2 "command: ${cmd}"
	eval "$cmd"
	rc=$?
	end=$(date +%s)
	local dur=$((end - start))
	if ((rc == 0)); then
		polap_log1 "✔ step: ${name} (took ${dur}s)"
	else
		polap_log0 "✖ step: ${name} FAILED (rc=${rc}, ${dur}s)"
		return "$rc"
	fi
}
