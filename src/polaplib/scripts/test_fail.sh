#!/usr/bin/env bash
# Version: v0.2.0

set -Eeuo pipefail
set -o errtrace
set -o functrace

POLAPLIB="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${POLAPLIB}/polap-lib-crash.sh"
crash_enable

message="intentional test failure (bash)"
exit_code=7
infile=""
outfile=""

# Robust long-option parser
while [[ $# -gt 0 ]]; do
	case "$1" in
	--message)
		message="${2:-$message}"
		shift 2
		;;
	--exit-code)
		exit_code="${2:-$exit_code}"
		shift 2
		;;
	-h | --help)
		printf 'Usage: bash %s [--message STR] [--exit-code N] [INFILE [OUTFILE]]\n' "${BASH_SOURCE[0]}"
		exit 0
		;;
	--)
		shift
		break
		;;
	-*)
		printf '[%s] ERROR: unknown option %s\n' "$(_here)" "$1" >&2
		exit 2
		;;
	*)
		if [[ -z "$infile" ]]; then infile="$1"; else outfile="$1"; fi
		shift
		;;
	esac
done

# Optional write
if [[ -n "$outfile" ]]; then
	printf 'test_fail.sh wrote before failing\n' >"$outfile"
fi

printf '[%s] scripts test bash: message=%s code=%s in=%s out=%s\n' \
	"$(_here)" "$message" "$exit_code" "${infile:-}" "${outfile:-}" >&2

third() {
	bash -c "exit ${exit_code}"
} # failing simple command

second() {
	third
}

first() {
	second
}

first
