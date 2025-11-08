#!/usr/bin/env bash
# polaplib/tooltest/minimap2/smoke.sh
set -Eeuo pipefail
tdir="$(cd "$(dirname "$0")" && pwd)"
# shellcheck disable=SC1091
source "${tdir}/../tooltest-lib.sh"
junit_begin

smoke() {
	# Check version in a safe window; relax as needed.
	assert_version_range "minimap2 --version" "2.26" "2.31"

	local ref="${tdir}/fixtures/ref.fa"
	local reads="${tdir}/fixtures/reads.fq"

	run_argv_to "minimap2/smoke.paf" -- minimap2 -x map-ont "$ref" "$reads"
	assert_file "${TOOLTEST_OUTDIR}/minimap2/smoke.paf"
}
run_case smoke smoke
junit_end
