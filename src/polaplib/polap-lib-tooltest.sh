#!/usr/bin/env bash
# polaplib/polap-lib-tooltest.sh
# Tiny helpers for command-line tool smoke tests.

# shellcheck disable=SC1091
source "${BASH_SOURCE[0]%/*}/polap-lib-run-common.sh"

: "${POLAP_TOOLTEST_OUTDIR:=o/tooltest}"

ptt_banner() { polap_log1 "=== $* ==="; }
ptt_require_tool() {
	local bin="$1"
	if ! command -v "$bin" >/dev/null 2>&1; then
		polap_log0 "SKIP: '$bin' not in PATH"
		return 1
	fi
	return 0
}
ptt_mkout() { polap__ensure_dir "$1"; }
ptt_assert_file() { polap_assert_file "$1"; }
ptt_assert_contains() { polap_assert_contains "$1" "$2"; }

# minimal fixtures
ptt_mini_fasta() {
	local f="$1"
	cat >"$f" <<'FA'
>refA
ACGTACGTACGTACGTACGT
>refB
TTTTAAAACCCCGGGGTTAA
FA
}
ptt_mini_fastq() {
	local f="$1"
	cat >"$f" <<'FQ'
@r1
ACGTACGTACGT
+
FFFFFFFFFFFF
@r2
TTTTAAAACCCC
+
FFFFFFFFFFFF
FQ
}
