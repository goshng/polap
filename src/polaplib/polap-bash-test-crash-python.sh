#!/usr/bin/env bash
# polaplib/polap-bash-test-crash-python.sh
# Version: v0.1.0
#
# PURPOSE
#   Minimal parent launcher to exercise a Python crash in:
#     polaplib/scripts/test_crash.py
#   This parent enables the Bash failsafe, then invokes the child so you get:
#     • the Python traceback from the child, and
#     • (if the non-zero bubbles up and failsafe is active here) a one-liner
#       from this parent showing the call site.
#
# USAGE
#   bash polaplib/polap-bash-test-crash-python.sh [--message STR] [--code N] [IN [OUT]]
#
# EXAMPLES
#   bash polaplib/polap-bash-test-crash-python.sh
#   bash polaplib/polap-bash-test-crash-python.sh --message "provoking" --code 9 in.txt out.txt

# --- ensure bash / safe re-exec under sh --------------------------------------
had_u=0
case $- in *u*) had_u=1 ;; esac
set +u
if [ -z "${BASH_VERSION+x}" ]; then exec /usr/bin/env bash "$0" "$@"; fi
[ "$had_u" -eq 1 ] && set -u

set -Eeuo pipefail
set -o errtrace
set -o functrace

# --- locate POLAPLIB and enable failsafe for THIS process ---------------------
POLAPLIB="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${POLAPLIB}/polap-lib-failsafe.sh" # your v0.3.0 baseline
polap_enable_failsafe

# --- defaults for arguments we forward to Python ------------------------------
message="intentional test failure (python)"
code=7
infile=""
outfile=""

# --- tiny argv parser (Bash) --------------------------------------------------
while [[ $# -gt 0 ]]; do
	case "$1" in
	--message)
		message="${2:-$message}"
		shift 2
		;;
	--code | --exit-code)
		code="${2:-$code}"
		shift 2
		;;
	-h | --help)
		cat <<EOF
Usage: $0 [--message STR] [--code N] [INFILE [OUTFILE]]
EOF
		exit 0
		;;
	--)
		shift
		break
		;;
	-*)
		echo "[ERROR] unknown option: $1" >&2
		exit 2
		;;
	*)
		if [[ -z "$infile" ]]; then infile="$1"; else outfile="$1"; fi
		shift
		;;
	esac
done

# --- child path ---------------------------------------------------------------
child_py="${POLAPLIB}/scripts/test_crash.py"

# --- context header from parent (optional) ------------------------------------
printf '[%s:%s] parent(py): message=%s code=%s in=%s out=%s\n' \
	"${BASH_SOURCE[0]##*/}" "${LINENO}" "$message" "$code" "${infile:-}" "${outfile:-}" >&2

# --- execute the child (this is where the crash should happen) ----------------
# If the child exits non-zero and set -e is active, our ERR trap (failsafe) will
# print this call site.
python3 "$child_py" --message "$message" --exit-code "$code" ${infile:+--in "$infile"} ${outfile:+--out "$outfile"}
