#!/usr/bin/env bash
# polaplib/polap-bash-tooltest-run-all.sh
#
# Run all tooltests (minimap2, samtools, iqtree2, flye), aggregate PASS/SKIP/FAIL,
# and print a summary table. Assumes each test lives at:
#   polaplib/tooltest/<tool>/run.sh
#
# Exit code: 0 if no FAIL; 1 if any FAIL; 99 if everything skipped.
#
# Options:
#   --list              List discovered tests and exit
#   --pattern REGEX     Only run tools whose name matches REGEX
#   --verbose N         Set POLAP_VERBOSE (0..4)
#   --outdir DIR        Override POLAP_TOOLTEST_OUTDIR for children
#
# Examples:
#   bash polaplib/polap-bash-tooltest-run-all.sh
#   bash polaplib/polap-bash-tooltest-run-all.sh --pattern 'mini|sam'
#   bash polaplib/polap-bash-tooltest-run-all.sh --list
#
set -Eeuo pipefail

# ── Locate polaplib and common logging helpers ────────────────────────────────
BASE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${BASE}/polap-lib-run-common.sh"

# Default log file if none provided
: "${LOG_FILE:=o/tooltest/run-all.log}"
polap__mkparent "$LOG_FILE"

# Defaults passed down to tests (they also use this by convention)
: "${POLAP_TOOLTEST_OUTDIR:=o/tooltest}"

# ── CLI parsing (simple) ─────────────────────────────────────────────────────
LIST_ONLY=0
PATTERN=""
while [[ $# -gt 0 ]]; do
	case "$1" in
	--list)
		LIST_ONLY=1
		shift
		;;
	--pattern)
		PATTERN="${2:-}"
		shift 2
		;;
	--verbose)
		export POLAP_VERBOSE="${2:-1}"
		shift 2
		;;
	--outdir)
		export POLAP_TOOLTEST_OUTDIR="${2:-$POLAP_TOOLTEST_OUTDIR}"
		shift 2
		;;
	--)
		shift
		break
		;;
	*)
		# Treat bare args as tool names to run (override default list)
		break
		;;
	esac
done

# Tools to run (default) or from remaining args
declare -a TOOLS_DEFAULT=(minimap2 samtools iqtree2 flye)
declare -a TOOLS=()
if [[ $# -gt 0 ]]; then
	TOOLS=("$@")
else
	TOOLS=("${TOOLS_DEFAULT[@]}")
fi

# ── Discover tests ───────────────────────────────────────────────────────────
declare -a NAMES=() PATHS=()
for t in "${TOOLS[@]}"; do
	p="${BASE}/tooltest/${t}/run.sh"
	if [[ -f "$p" ]]; then
		if [[ -n "$PATTERN" ]]; then
			if [[ "$t" =~ $PATTERN ]]; then
				NAMES+=("$t")
				PATHS+=("$p")
			fi
		else
			NAMES+=("$t")
			PATHS+=("$p")
		fi
	else
		polap_log0 "MISSING: ${t} (no ${p})"
	fi
done

if ((LIST_ONLY)); then
	if ((${#NAMES[@]} == 0)); then
		polap_log0 "No tests discovered."
		exit 0
	fi
	polap_log1 "Discovered tests:"
	for i in "${!NAMES[@]}"; do
		printf '  - %-10s  %s\n' "${NAMES[$i]}" "${PATHS[$i]}"
	done
	exit 0
fi

if ((${#NAMES[@]} == 0)); then
	polap_log0 "Nothing to run."
	exit 0
fi

# ── Runner (does not abort on child failure) ─────────────────────────────────
_status_of_rc() {
	local rc="$1"
	if ((rc == 0)); then
		printf 'PASS'
	elif ((rc == 99)); then
		printf 'SKIP'
	else
		printf 'FAIL'
	fi
}

# Arrays for report
declare -a STATS=() RCS=() DURS=() FULLS=()

TOTAL_PASS=0
TOTAL_SKIP=0
TOTAL_FAIL=0

polap_log1 "Running ${#NAMES[@]} tool test(s) → outdir=${POLAP_TOOLTEST_OUTDIR}"

for i in "${!NAMES[@]}"; do
	name="${NAMES[$i]}"
	path="${PATHS[$i]}"

	# Preview
	polap_log2 "command: bash $(printf %q "$path")"

	# Time + capture rc without tripping set -e
	start=$(date +%s)
	set +e
	bash "$path"
	rc=$?
	set -e
	end=$(date +%s)
	dur=$((end - start))

	status="$(_status_of_rc "$rc")"
	case "$status" in
	PASS) ((TOTAL_PASS++)) ;;
	SKIP) ((TOTAL_SKIP++)) ;;
	FAIL) ((TOTAL_FAIL++)) ;;
	esac

	STATS+=("$status")
	RCS+=("$rc")
	DURS+=("$dur")
	FULLS+=("$path")

	polap_log1 "$(printf '→ %-10s  %-4s (rc=%d, %ds)' "$name" "$status" "$rc" "$dur")"
done

# ── Summary table ────────────────────────────────────────────────────────────
echo ""
printf '%s\n' "Summary:"
printf '%-12s %-6s %4s %6s  %s\n' "tool" "status" "rc" "time" "script"
printf '%-12s %-6s %4s %6s  %s\n' "------------" "------" "----" "------" "-------------------------------"
for i in "${!NAMES[@]}"; do
	printf '%-12s %-6s %4d %6ss  %s\n' \
		"${NAMES[$i]}" "${STATS[$i]}" "${RCS[$i]}" "${DURS[$i]}" "${FULLS[$i]}"
done
echo ""
printf 'PASS: %d   SKIP: %d   FAIL: %d\n' "$TOTAL_PASS" "$TOTAL_SKIP" "$TOTAL_FAIL"

# Exit code policy
if ((TOTAL_FAIL > 0)); then
	exit 1
elif ((TOTAL_PASS == 0)) && ((TOTAL_SKIP > 0)); then
	# everything skipped
	exit 99
else
	exit 0
fi
