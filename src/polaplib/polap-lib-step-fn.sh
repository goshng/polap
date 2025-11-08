#!/usr/bin/env bash
# polaplib/polap-lib-step-fn.sh
# polap_step_fn: run a bash function as a named step.

# shellcheck disable=SC1091
source "${BASH_SOURCE[0]%/*}/polap-lib-run-common.sh"

# Usage: polap_step_fn "<step-name>" function_name [args...]
polap_step_fn() {
	local name="$1" fn="$2"
	shift 2
	local start end rc
	start=$(date +%s)
	polap_log1 "▶ step(fn): ${name} → ${fn} $*"
	"$fn" "$@"
	rc=$?
	end=$(date +%s)
	local dur=$((end - start))
	if ((rc == 0)); then
		polap_log1 "✔ step(fn): ${name} (took ${dur}s)"
	else
		polap_log0 "✖ step(fn): ${name} FAILED (rc=${rc}, ${dur}s)"
		return "$rc"
	fi
}
