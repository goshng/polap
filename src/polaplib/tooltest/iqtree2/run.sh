#!/usr/bin/env bash
# polaplib/tooltest/iqtree2/run.sh
# Tiny alignment → infer tree. Uses iqtree2 if present; falls back to iqtree.

set -Eeuo pipefail

TT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${TT_DIR}/../common.sh"

OUT="${POLAP_TOOLTEST_OUTDIR}/iqtree2"
tt_mkout "${OUT}"
tt_banner "iqtree2 tiny alignment test"

# Resolve binary: prefer iqtree2, fallback to iqtree
IQBIN=""
if command -v iqtree2 >/dev/null 2>&1; then
	IQBIN="iqtree2"
elif command -v iqtree >/dev/null 2>&1; then
	IQBIN="iqtree"
else
	printf '[%s %s%s] SKIP: iqtree2/iqtree not found in PATH\n' "$(tt_ts)" "$(tt_conda_tag)" "${FUNCNAME[0]}" >&2
	exit 99
fi

# Fixture: tiny DNA alignment
ALN="${OUT}/aln.fa"
tt_ensure_mini_aln_fa "${ALN}"

# Clean any previous run with the same prefix
rm -f "${OUT}/demo."*

# Run (small, deterministic)
# -nt 1 to be safe on tiny inputs; -seed 1 for reproducibility; -quiet to reduce noise
tt_run "${IQBIN}" -s "${ALN}" -m GTR+G -nt 1 -seed 1 -pre "${OUT}/demo" -quiet

# Check outputs
TREE="${OUT}/demo.treefile"
LOG="${OUT}/demo.log"
IQREPORT="${OUT}/demo.iqtree"

tt_assert_file "${TREE}"
tt_assert_file "${IQREPORT}"
tt_assert_contains "(" "${TREE}"

tt_banner "iqtree2 OK → ${OUT}"
