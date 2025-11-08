#!/usr/bin/env bash
# polaplib/scripts/test_crash.sh
# Version: v0.1.1
# Intentionally fails in different ways so you can verify the failsafe.

# ensure bash, keep nounset safe across re-exec
had_u=0
case $- in *u*) had_u=1 ;; esac
set +u
if [ -z "${BASH_VERSION+x}" ]; then exec /usr/bin/env bash "$0" "$@"; fi
[ "$had_u" -eq 1 ] && set -u

set -Eeuo pipefail
set -o errtrace
set -o functrace

POLAPLIB="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck disable=SC1091
source "${POLAPLIB}/polap-lib-failsafe.sh"
polap_enable_failsafe

mode="cmd"
code=7
while [[ $# -gt 0 ]]; do
	case "$1" in
	--mode)
		mode="${2:-$mode}"
		shift 2
		;;
	--code)
		code="${2:-$code}"
		shift 2
		;;
	-h | --help)
		echo "Usage: $0 [--mode cmd|exit|nounset] [--code N]"
		exit 0
		;;
	*) break ;;
	esac
done

test_crash_fn1() {
	bash -c "exit ${code}"
}

printf '[%s:%s] test_crash: mode=%s code=%s\n' "${BASH_SOURCE[0]##*/}" "${LINENO}" "$mode" "$code" >&2

case "$mode" in
cmd)
	test_crash_fn1
	# bash -c "exit ${code}"
	;; # preferred for ERR
exit)
	exit "${code}"
	;; # parent’s failsafe will report the call-site
nounset)
	unset __X
	echo "$__X"
	;; # trigger set -u
*)
	echo "[ERROR] unknown --mode '$mode'" >&2
	exit 2
	;;
esac
