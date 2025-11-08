#!/usr/bin/env bash
# polaplib/polap-bash-test-crash-bash.sh
# Version: v0.1.1
#
# PURPOSE
#   Thin parent launcher that:
#     1) turns on the failsafe for *itself*,
#     2) invokes the child crash script with bash,
#     3) lets both parent and child report their own sites if they fail.
#
# USAGE
#   bash polaplib/polap-bash-test-crash-bash.sh [--mode ...] [--code N]
#
# EXAMPLES
#   # failing simple command in the child:
#   bash polaplib/polap-bash-test-crash-bash.sh --mode cmd --code 9
#
#   # nounset in the child:
#   bash polaplib/polap-bash-test-crash-bash.sh --mode nounset
#
#   # direct exit in the child:
#   bash polaplib/polap-bash-test-crash-bash.sh --mode exit --code 5

had_u=0
case $- in *u*) had_u=1 ;; esac
set +u
if [ -z "${BASH_VERSION+x}" ]; then exec /usr/bin/env bash "$0" "$@"; fi
[ "$had_u" -eq 1 ] && set -u

set -Eeuo pipefail
set -o errtrace
set -o functrace

POLAPLIB="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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

child="${POLAPLIB}/scripts/test_crash.sh"
printf '[%s:%s] parent: mode=%s code=%s\n' "${BASH_SOURCE[0]##*/}" "${LINENO}" "$mode" "$code" >&2

# child’s own failsafe will report failures at its site;
# this parent’s failsafe will report the call-site here if the child exits non-zero.
# bash "$child" --mode "$mode" --code "$code"
bash "$child" --mode "$mode" --code "$code"
