#!/usr/bin/env bash
# polaplib/polap-lib-crash.sh
# Version: v0.1.2 — nounset-safe (_here, trap, stack)

set -o errtrace
set -o functrace

# ----- where-am-I helper (func:file:line), nounset-safe -----
_here() {
	if caller 0 >/dev/null 2>&1; then
		# caller prints: "<line> <func> <file>"
		if read -r __l __f __fl < <(caller 0); then
			: "${__f:=main}"
			printf '%s:%s:%s' "$__f" "${__fl##*/}" "$__l"
			return
		fi
	fi
	printf '%s:%s:%s' main "${BASH_SOURCE[0]##*/}" "${LINENO:-0}"
}

# ----- compact stack printer using `caller`, nounset-safe -----
print_call_chain() {
	# NEVER rely on an uninitialized i under set -u
	i=0
	ts="$(date '+%Y-%m-%d %H:%M:%S')"
	# use ${i:-0} everywhere we read it
	while caller "${i:-0}" >/dev/null 2>&1; do
		if read -r l f fl < <(caller "${i:-0}"); then
			: "${f:=main}"
			src=""
			if [[ -r "${fl:-}" && "${l:-}" =~ ^[0-9]+$ ]]; then
				src="$(sed -n "${l}p" -- "$fl" 2>/dev/null | tr -d $'\r' | tr -s '[:space:]' ' ')"
			fi
			if [[ -n "$src" ]]; then
				printf '[%s %s@%s:%s] CALLED FROM: %s\n' \
					"$ts" "$f" "${fl##*/}" "$l" "$src" >&2
			else
				printf '[%s %s@%s:%s] CALLED FROM\n' \
					"$ts" "$f" "${fl##*/}" "$l" >&2
			fi
		fi
		i=$((${i:-0} + 1))
	done
}

# ----- lightweight crash enable (ERR only), nounset-safe -----
crash_enable() {

	trap '
    __st=$?; set +u
    __cmd=${BASH_COMMAND-}
    # prefer caller 0; fall back to arrays
    if read -r __l __f __fl < <(caller 0); then
      : "${__f:=main}"
      __file="$__fl"; __line="$__l"; __func="$__f"
    else
      __file="${BASH_SOURCE[1]-${BASH_SOURCE[0]-$0}}"
      __line="${BASH_LINENO[0]-${LINENO-0}}"
      __func="${FUNCNAME[1]-main}"
    fi
    printf "[%(%Y-%m-%d %H:%M:%S)T %s@%s:%s] FAIL: %s (exit %d)\n" \
           -1 "$__func" "${__file##*/}" "$__line" \
           "${__cmd//$'\n'/\\n}" "$__st" >&2
    print_call_chain
    set -u
  ' ERR

	# Optional per-command trace (toggle with POLAP_TRACE=1)
	if [[ "${POLAP_TRACE:-0}" -eq 1 ]]; then
		export BASH_XTRACEFD=9
		exec 9> >(awk '{ t=strftime("[%Y-%m-%d %H:%M:%S]"); print t, $0 }' >&2)
		export PS4='[${FUNCNAME[0]-main}@${BASH_SOURCE##*/}:${LINENO}] '
		set -x
	fi
}
