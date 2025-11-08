#!/usr/bin/env bash
# polaplib/polap-lib-step-idem.sh
# polap_step_idem: skip if marker exists; otherwise run string cmd and create marker.

# shellcheck disable=SC1091
source "${BASH_SOURCE[0]%/*}/polap-lib-run-common.sh"

# Usage: polap_step_idem "<step-name>" "<marker-file>" "<cmd string>"
polap_step_idem() {
	local name="$1" marker="$2" cmd="$3"
	if [[ -e "$marker" ]]; then
		polap_log1 "↺ step(idem): ${name} (skip; marker=${marker})"
		return 0
	fi
	polap__mkparent "$marker"
	polap_log1 "▶ step(idem): ${name}"
	polap_log2 "command: ${cmd}"
	local start end rc
	start=$(date +%s)
	eval "$cmd"
	rc=$?
	end=$(date +%s)
	local dur=$((end - start))
	if ((rc == 0)); then
		: >"$marker"
		polap_log1 "✔ step(idem): ${name} (took ${dur}s)"
	else
		polap_log0 "✖ step(idem): ${name} FAILED (rc=${rc}, ${dur}s)"
		return "$rc"
	fi
}
