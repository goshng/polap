#!/usr/bin/env bash
# polaplib/tooltest/minimap2/gold.sh
# This "gold" test is semantic (stable fields) rather than byte-equal.
set -Eeuo pipefail
tdir="$(cd "$(dirname "$0")" && pwd)"
# shellcheck disable=SC1091
source "${tdir}/../tooltest-lib.sh"
junit_begin

gold_semantic() {
	local ref="${tdir}/fixtures/ref.fa"
	local reads="${tdir}/fixtures/reads.fq"
	local out="${TOOLTEST_OUTDIR}/minimap2/gold.paf"

	run_argv_to "minimap2/gold.paf" -- minimap2 -x map-ont "$ref" "$reads"
	assert_file "$out"

	# Expect exactly 1 alignment (our fixture is trivial).
	local n
	n="$(grep -vc '^#' "$out" || true)"
	[[ "$n" -eq 1 ]] || {
		logi "FAIL: expected 1 PAF line, got $n"
		return 1
	}

	# Check key fields: qname, strand, tname, qlen consistency
	awk -F'\t' '
    BEGIN{ok=1}
    $0 !~ /^#/ {
      q=$1; ql=$2+0; qs=$3+0; qe=$4+0; st=$5; t=$6;
      if(q!="read1"){print "bad qname: " q > "/dev/stderr"; ok=0}
      if(st!="+" && st!="-"){print "bad strand: " st > "/dev/stderr"; ok=0}
      if(t!="chr1"){print "bad tname: " t > "/dev/stderr"; ok=0}
      if(qe-qs!=ql){print "inconsistent q span vs qlen: " qe-qs " vs " ql > "/dev/stderr"; ok=0}
    }
    END{exit ok?0:1}
  ' "$out"
}
run_case gold_semantic gold_semantic
junit_end
