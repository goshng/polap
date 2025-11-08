#!/usr/bin/env bash
# polaplib/polap-lib-step-argv.sh
# polap_step_argv: named step wrapper around argv, with timing and logs.

# shellcheck disable=SC1091
source "${BASH_SOURCE[0]%/*}/polap-lib-run-common.sh"

# Usage: polap_step_argv "<step-name>" cmd arg1 arg2 ...
polap_step_argv() {
	local name="$1"
	shift
	local start end rc
	start=$(date +%s)

	polap_log1 "▶ step: ${name}"
	local preview
	preview="$(polap__preview_argv "$@")"
	polap_log2 "command: ${preview}"

	"$@"
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
