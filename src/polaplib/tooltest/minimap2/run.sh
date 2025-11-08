#!/usr/bin/env bash
# polaplib/tooltest/minimap2/run.sh

set -Eeuo pipefail

TT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="${TT_DIR%/tooltest/*}" # optional, if you want a base
# shellcheck disable=SC1091
source "${ROOT_DIR}/polaplib/polap-lib-tooltest.sh"
# shellcheck disable=SC1091
source "${ROOT_DIR}/polaplib/polap-lib-run-argv.sh"

OUT="${POLAP_TOOLTEST_OUTDIR}/minimap2"
ptt_mkout "${OUT}"

ptt_banner "minimap2 tiny mapping"

if ! ptt_require_tool minimap2; then
	exit 99
fi

REF="${OUT}/ref.fa"
READS="${OUT}/reads.fq"
SAM="${OUT}/map.sam"

ptt_mini_fasta "${REF}"
ptt_mini_fastq "${READS}"

# Use a simple preset; FASTA ref + FASTQ reads is fine with general map preset.
polap_run_argv minimap2 -a "${REF}" "${READS}" >"${SAM}"

ptt_assert_file "${SAM}"
ptt_assert_contains "@SQ" "${SAM}" # header present
ptt_assert_contains "@PG" "${SAM}" # program line
# At least one alignment line (non-header)
if ! grep -v '^@' "${SAM}" | head -n1 | grep -q .; then
	polap_log0 "No alignments produced in ${SAM}"
	exit 1
fi

ptt_banner "minimap2 OK → ${OUT}"
