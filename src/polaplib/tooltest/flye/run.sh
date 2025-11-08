#!/usr/bin/env bash
# polaplib/tooltest/flye/run.sh
# Smoke-only test: verify binary and show version/help.
# (Avoids heavy assembly; you can enable a heavier test via env if needed.)

set -Eeuo pipefail

TT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${TT_DIR}/../_lib/common.sh"

OUT="${POLAP_TOOLTEST_OUTDIR}/flye"
tt_mkout "${OUT}"

tt_banner "flye smoke test (version/help)"

# Require tool
tt_require_tool flye || exit 99

# Version
VERS_TXT="${OUT}/version.txt"
tt_run_to "${VERS_TXT}" "${OUT}/version.err" flye --version
tt_assert_file "${VERS_TXT}"

# Help (first lines)
HELP_TXT="${OUT}/help.txt"
# On some versions, `--help` prints to stdout; capture either way
tt_run bash -c 'flye --help 2>&1 | head -n 30' >"${HELP_TXT}"
tt_assert_file "${HELP_TXT}"

tt_banner "flye OK (smoke) → ${OUT}"

# Optional heavier check (OFF by default): set POLAP_TOOLTEST_FLYE_RUN=1 to try
# a tiny run. NOTE: This is likely to fail on the miniature fixtures; it is
# primarily useful to confirm error handling and logging end-to-end.
if [[ "${POLAP_TOOLTEST_FLYE_RUN:-0}" == "1" ]]; then
	tt_banner "flye optional tiny run (may fail by design)"
	READS="${OUT}/mini.fq"
	tt_ensure_mini_fq "${READS}"
	ASM="${OUT}/asm"
	rm -rf -- "${ASM}"
	# This may exit non-zero; do NOT let it crash the suite.
	set +e
	tt_run flye --nano-raw "${READS}" --genome-size 100k --threads 1 --out-dir "${ASM}"
	rc=$?
	set -e
	printf '[%s %s%s] flye tiny run exited rc=%d (expected may be nonzero)\n' "$(tt_ts)" "$(tt_conda_tag)" "${FUNCNAME[0]}" "$rc" >&2
fi
