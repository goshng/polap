#!/usr/bin/env bash
# polaplib/polap-lib-run.sh
# polap_run: run a *string* command via eval, log literal + expanded forms.

# shellcheck disable=SC1091
source "${BASH_SOURCE[0]%/*}/polap-lib-run-common.sh"

polap_run_simple() {
	# Usage: polap_run_simple [--tag TAG] -- "<single shell command string>"
	local tag="" cmd=""
	while [[ $# -gt 0 ]]; do
		case "$1" in
		--tag)
			tag="$2"
			shift 2
			;;
		--)
			shift
			break
			;;
		*) break ;;
		esac
	done
	# The rest is the single command string (joined with spaces)
	cmd="$*"
	_polap_log0 "[RUN:${tag}] ${cmd}"
	# Optional: show parameter-expanded form (no command substitution), if you have it
	# _polap_log0 "[EXPAND:${tag}] $(polap__expand_string "$cmd")"
	eval "$cmd"
}
