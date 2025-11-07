#!/usr/bin/env bash
# Version: v0.2.0

set -Eeuo pipefail
set -o errtrace
set -o functrace

POLAPLIB="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck disable=SC1091
source "${POLAPLIB}/polap-lib-crash.sh"
crash_enable

BASH_FAIL="${POLAPLIB}/scripts/test_fail.sh"
PY_FAIL="${POLAPLIB}/scripts/test_fail.py"
R_FAIL="${POLAPLIB}/scripts/test_fail.R"

which="bash"
exit_code=7
message="provoking nested failures across runtimes"
infile="infile.txt"
outfile="outfile.txt"

while [[ $# -gt 0 ]]; do
	case "$1" in
	--which)
		which="${2:-$which}"
		shift 2
		;;
	--exit-code)
		exit_code="${2:-$exit_code}"
		shift 2
		;;
	--message)
		message="${2:-$message}"
		shift 2
		;;
	-h | --help)
		echo "Usage: bash polap-bash-test-fail.sh [--which bash|python|r|all] [--exit-code N] [--message STR] [IN [OUT]]"
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
	*) break ;;
	esac
done
[[ $# -ge 1 ]] && infile="$1"
[[ $# -ge 2 ]] && outfile="$2"

printf '[%s] test: which=%s exit=%s message=%s in=%s out=%s\n' \
	"$(_here)" "$which" "$exit_code" "$message" "$infile" "$outfile" >&2

case "$which" in
bash)
	bash "$BASH_FAIL" --message "$message" --exit-code "$exit_code" "$infile" "$outfile"
	;;
python)
	python3 "$PY_FAIL" --message "$message" --exit-code "$exit_code" "$infile" "$outfile"
	;;
r)
	Rscript --vanilla "$R_FAIL" --message "$message" --exit-code "$exit_code" "$infile" "$outfile"
	;;
all)
	bash "$BASH_FAIL" --message "$message" --exit-code "$exit_code" "$infile" "$outfile"
	python3 "$PY_FAIL" --message "$message" --exit-code "$exit_code" "$infile" "$outfile"
	Rscript --vanilla "$R_FAIL" --message "$message" --exit-code "$exit_code" "$infile" "$outfile"
	;;
*)
	echo "[ERROR] --which must be bash|python|r|all" >&2
	exit 2
	;;
esac
