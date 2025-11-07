#!/usr/bin/env bash
# Version: v0.2.0

set -Eeuo pipefail
set -o errtrace
set -o functrace

POLAPLIB="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${POLAPLIB}/polap-lib-crash.sh"
crash_enable

cmd_test_main() {
	local which="${1:-bash}" code="${2:-7}" msg="${3:-provoking nested failures across runtimes}"
	local infile="${4:-infile.txt}" outfile="${5:-outfile.txt}"

	printf '[%s] Running test: which=%s exit=%s message=%s in=%s out=%s\n' \
		"$(_here)" "$which" "$code" "$msg" "$infile" "$outfile" >&2

	bash "${POLAPLIB}/polap-bash-test-fail.sh" \
		--which "$which" --exit-code "$code" --message "$msg" "$infile" "$outfile"
}

cmd_test_py() {
	local code=7
	local msg="message run python3"
	local infile="infile.txt"
	local outfile="outfile.txt"

	python3 "${POLAPLIB}/polap-py-test-fail.py" --exit-code "$code" --message "$msg" "$infile" "$outfile"
}

cmd_test_r() {
	local code=7
	local msg="message run R"
	local infile="infile.txt"
	local outfile="outfile.txt"
	Rscript --vanilla "${POLAPLIB}/polap-r-test-fail.R" --exit-code "$code" --message "$msg" "$infile" "$outfile"

}
