#!/usr/bin/env bash
set -Eeuo pipefail

TT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
POLAPLIB="$(cd "${TT_DIR}/../.." && pwd)"
# shellcheck disable=SC1091
source "${POLAPLIB}/tooltest/_lib/common.sh"

OUT="${POLAP_TOOLTEST_OUTDIR}/flye"
tt_mkout "${OUT}"

EXTENDED=0
while [[ $# -gt 0 ]]; do
	case "$1" in
	--extended)
		EXTENDED=1
		shift
		;;
	-h | --help)
		echo "Usage: $(basename "$0") [--extended]"
		exit 0
		;;
	*)
		echo "Unknown arg: $1" >&2
		exit 2
		;;
	esac
done

tt_banner "flye: smoke (version/help)"
tt_require_tool flye || exit 99

VER="${OUT}/version.txt"
tt_run_to "$VER" "${OUT}/version.err" flye --version
tt_assert_file "$VER"

HELP="${OUT}/help.txt"
# Help sometimes goes to stdout; normalize
tt_run bash -c 'flye --help 2>&1 | head -n 30' >"$HELP"
tt_assert_file "$HELP"

if ((EXTENDED)); then
	tt_banner "flye: extended tiny run (may fail by design on micro data)"
	READS="${OUT}/mini.fq"
	ASM="${OUT}/asm"
	tt_ensure_mini_fq "$READS"
	rm -rf -- "$ASM"
	set +e
	tt_run flye --nano-raw "$READS" --genome-size 100k --threads 1 --out-dir "$ASM"
	rc=$?
	set -e
	printf '[%s %s] flye tiny run exited rc=%d (nonzero is acceptable in extended)\n' "$(tt_ts)" "$(tt_conda_tag)" "$rc" >&2
fi

tt_banner "flye: OK → ${OUT}"
