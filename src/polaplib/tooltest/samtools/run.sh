#!/usr/bin/env bash
# polaplib/tooltest/samtools/run.sh
# Smoke + tiny workflow test for samtools using a minimal SAM fixture.

set -Eeuo pipefail

# Source shared helpers
TT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${TT_DIR}/../common.sh"

OUT="${POLAP_TOOLTEST_OUTDIR}/samtools"
tt_mkout "${OUT}"

tt_banner "samtools smoke + convert/sort/index"

# Require tool
tt_require_tool samtools || exit 99

# Fixtures: minimal SAM
SAM="${OUT}/mini.sam"
tt_ensure_mini_sam "${SAM}"

# 1) Convert SAM → BAM
RAW_BAM="${OUT}/mini.bam"
tt_run samtools view -b -o "${RAW_BAM}" "${SAM}"
tt_assert_file "${RAW_BAM}"

# 2) Sort BAM
SORTED_BAM="${OUT}/mini.sorted.bam"
tt_run samtools sort -o "${SORTED_BAM}" "${RAW_BAM}"
tt_assert_file "${SORTED_BAM}"

# 3) Index BAM
tt_run samtools index "${SORTED_BAM}"
tt_assert_file "${SORTED_BAM}.bai"

# 4) Quick stats
FLAGSTAT_TXT="${OUT}/flagstat.txt"
IDXSTATS_TXT="${OUT}/idxstats.txt"
tt_run_to "${FLAGSTAT_TXT}" "${OUT}/flagstat.err" samtools flagstat "${SORTED_BAM}"
tt_run_to "${IDXSTATS_TXT}" "${OUT}/idxstats.err" samtools idxstats "${SORTED_BAM}"

tt_assert_file "${FLAGSTAT_TXT}"
tt_assert_file "${IDXSTATS_TXT}"
tt_assert_contains "in total" "${FLAGSTAT_TXT}"

tt_banner "samtools OK → ${OUT}"
