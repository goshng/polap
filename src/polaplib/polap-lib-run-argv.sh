#!/usr/bin/env bash
# polaplib/polap-lib-run-argv.sh
# polap_run_argv: run argv safely, no eval.

# shellcheck disable=SC1091
source "${BASH_SOURCE[0]%/*}/polap-lib-run-common.sh"

polap_run_argv() {
	# Usage: polap_run_argv [--tag TAG] -- cmd arg1 arg2 ...
	local tag=""
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
	local preview
	preview="$(printf '%q ' "$@")"
	preview="${preview# }"
	_polap_log0 "[RUN:${tag}] ${preview}"
	"$@"
}
