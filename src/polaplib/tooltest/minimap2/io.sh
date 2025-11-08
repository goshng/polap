#!/usr/bin/env bash
# polaplib/tooltest/minimap2/io.sh
set -Eeuo pipefail
tdir="$(cd "$(dirname "$0")" && pwd)"
# shellcheck disable=SC1091
source "${tdir}/../tooltest-lib.sh"
junit_begin

io_contract() {
	local ref="${tdir}/fixtures/ref.fa"
	local reads="${tdir}/fixtures/reads.fq"
	local out="${TOOLTEST_OUTDIR}/minimap2/io.paf"

	run_argv_to "minimap2/io.paf" -- minimap2 -x map-ont "$ref" "$reads"
	assert_file "$out"

	# PAF must have >=12 fields; check all non-comment lines
	awk 'BEGIN{bad=0} $0 !~ /^#/ && NF<12{bad=1} END{exit bad}' "$out"
}
run_case io_contract io_contract
junit_end
