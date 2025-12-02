#!/usr/bin/env bash
# File: polap-bash-run-command-selftest.sh
# Version: v0.2.4
# Purpose: End-to-end sanity checks for polap-lib-run-command.sh.

set -Eeuo pipefail

# Locate and source the runner
LIBDIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd -P)"
RUNLIB="${LIBDIR}/polap-lib-run-command.sh"
[[ -s "$RUNLIB" ]] || {
	echo "Missing: $RUNLIB" >&2
	exit 2
}
# shellcheck source=/dev/null
source "$RUNLIB"

# Fresh workspace
TMP="$(mktemp -d)"
export OUTDIR="${TMP}/out"
polap__ensure_dir "$OUTDIR"

passes=0 fails=0
pass() {
	echo "[PASS] $*"
	passes=$((passes + 1))
}
fail() {
	echo "[FAIL] $*" >&2
	fails=$((fails + 1))
}

echo "[selftest] OUTDIR=$OUTDIR"

# 1) argv backend (preferred; no eval) — /bin/true
POLAP_RUN_BACKEND=argv polap_run_wrapper "argv_ok" /bin/true
[[ -e "$OUTDIR/.polap-run/argv_ok/.ok" ]] && pass "argv_ok wrote .ok" || fail "argv_ok missing .ok"

# 2) argv backend skip works (idempotency)
POLAP_RUN_BACKEND=argv polap_run_wrapper "argv_ok" /bin/true && pass "argv_ok skipped when .ok exists"

# 3) failure propagation — /bin/false
if POLAP_RUN_BACKEND=argv polap_run_wrapper "argv_fail" /bin/false; then
	fail "argv_fail returned 0 unexpectedly"
else
	[[ ! -e "$OUTDIR/.polap-run/argv_fail/.ok" ]] && pass "argv_fail wrote no .ok on failure" || fail ".ok stamped on failure"
fi

POLAP_RUN_BACKEND=simple polap_run_wrapper "simple_pipe" "printf 'x\n' | wc -l"
[[ -e "$OUTDIR/.polap-run/simple_pipe/.ok" ]] && pass "simple_pipe wrote .ok" || fail "simple_pipe missing .ok"

# 5) direct backend
POLAP_RUN_BACKEND=direct polap_run_wrapper "direct_ok" echo "hello"
[[ -e "$OUTDIR/.polap-run/direct_ok/.ok" ]] && pass "direct_ok wrote .ok" || fail "direct_ok missing .ok"

# 6) dry-run mode (should not create .ok)
POLAP_DRYRUN=1 POLAP_RUN_BACKEND=argv polap_run_wrapper "dryrun_ok" /bin/true
[[ ! -e "$OUTDIR/.polap-run/dryrun_ok/.ok" ]] && pass "dryrun_ok did not write .ok" || fail "dryrun_ok incorrectly wrote .ok"
POLAP_DRYRUN=0

# 7) cmd.txt captures exact argv
POLAP_RUN_BACKEND=argv polap_run_wrapper "cmdtxt" printf "%s" "A B" "C"
if [[ -s "$OUTDIR/.polap-run/cmdtxt/cmd.txt" ]] && grep -q 'printf' "$OUTDIR/.polap-run/cmdtxt/cmd.txt"; then
	pass "cmd.txt recorded argv"
else
	fail "cmd.txt missing or incomplete"
fi

# 8) step.log captures output
POLAP_RUN_BACKEND=argv polap_run_wrapper "log_out" bash -c 'echo o1; echo e1 1>&2'
if [[ -s "$OUTDIR/.polap-run/log_out/step.log" ]] && grep -q 'o1' "$OUTDIR/.polap-run/log_out/step.log" && grep -q 'e1' "$OUTDIR/.polap-run/log_out/step.log"; then
	pass "step.log captured stdout/stderr"
else
	fail "step.log missing expected lines"
fi

# 9) lock helpers — acquire, contend, release
LOCKDIR="${OUTDIR}/.lock"
if polap_lock_acquire_wrapper "$LOCKDIR"; then
	pass "lock acquired first time"
else
	fail "lock acquisition failed unexpectedly"
fi
if polap_lock_acquire_wrapper "$LOCKDIR"; then
	fail "lock acquired twice (should have been busy)"
else
	pass "lock contention correctly reported non-zero"
fi
if polap_lock_release_wrapper "$LOCKDIR"; then
	pass "lock released ok"
else
	fail "lock release failed"
fi
# release on missing lock should still be 0
if polap_lock_release_wrapper "$LOCKDIR"; then
	pass "lock release on missing dir returns 0"
else
	fail "lock release on missing dir returned non-zero"
fi

# 10) logger never breaks pipelines (should not trigger 'else')
if polap_log1 "logger test"; then
	pass "logger returned 0"
else
	fail "logger returned non-zero"
fi

echo "[selftest] passes=$passes fails=$fails"
((fails == 0)) || exit 1
echo "[selftest] ALL PASS"
