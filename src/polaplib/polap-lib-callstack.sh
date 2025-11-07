#!/usr/bin/env bash
# Version: v0.3.0
# Minimal call-stack reporter that emits BOTH ERR and EXIT blocks,
# printing child + parent stacks, while de-duplicating identical prints.

# Config
: "${CS_RING_MAX:=64}"
: "${CS_MAX_WRAPPERS:=8}"
: "${CS_WRAPPER_REGEX:=^(bash|sh|python3?|Rscript)([[:space:]]|$)}"
# skip control/assign noise when choosing the “failed cmd”
: "${CS_SKIP_REGEX:=^([[:space:]]*(
    case(\s|$)|esac$|then$|fi$|do$|done$|elif(\s|$)|else$|
    \{|\}|;;|;|\(|\)|select(\s|$)|until(\s|$)|while(\s|$)|for(\s|$)|in(\s|$)|
    function(\s|$)|local(\s)|declare(\s)
  ).*|[[:space:]]*[A-Za-z_][A-Za-z0-9_]*=) }"

# ROLE/PID tag (export these at your call sites if you want custom labels)
: "${CS_ROLE:=PARENT}"
: "${CS_PID:=$$}"

# Ring
declare -a CS_CMD=() CS_FILE=() CS_LINE=() CS_FUNC=()
CS_IDX=0

cs_now() { date "+%Y-%m-%d %H:%M:%S"; }
cs_is_wrapper() { [[ "$1" =~ ${CS_WRAPPER_REGEX} ]]; }
cs_is_skippable() { [[ "$1" =~ ${CS_SKIP_REGEX} ]] || [[ -z "$1" ]]; }

cs_on_debug() {
	local cmd="${BASH_COMMAND-}"
	local line func file
	if read -r line func file < <(caller 0); then
		: "${func:=main}"
	else
		file="${BASH_SOURCE[1]-${BASH_SOURCE[0]}}"
		line="${BASH_LINENO[0]-${LINENO-0}}"
		func="${FUNCNAME[1]-main}"
	fi

	local i="$CS_IDX"
	CS_CMD[$i]="$cmd"
	CS_FILE[$i]="$file"
	CS_LINE[$i]="$line"
	CS_FUNC[$i]="$func"
	CS_IDX=$(((i + 1) % CS_RING_MAX))
}

cs_find_relevant_idx() {
	local head=$(((CS_IDX + CS_RING_MAX - 1) % CS_RING_MAX))
	local steps=0 i="$head"
	while ((steps < CS_RING_MAX)); do
		local c="${CS_CMD[$i]-}"
		if cs_is_skippable "$c"; then
			i=$(((i + CS_RING_MAX - 1) % CS_RING_MAX))
			((steps++))
			continue
		fi
		printf "%s" "$i"
		return 0
	done
	printf "%s" "$head"
}

# Dedupe identical print blocks
CS_LAST_SIG=""
cs_maybe_print_block() {
	# $1 = event (ERR|EXIT), $2 = file, $3 = line, $4 = cmd, $5 = ec
	local ev="$1" f="$2" l="$3" c="$4" ec="$5"
	local sig="${CS_ROLE}|${CS_PID}|${ev}|${f}|${l}|${c}|${ec}"
	if [[ "$sig" == "$CS_LAST_SIG" ]]; then
		return 1
	fi
	CS_LAST_SIG="$sig"
	return 0
}

cs_print_failed_at() {
	local ev="$1" idx="$2" ec="$3"
	local cmd="${CS_CMD[$idx]-}"
	local file="${CS_FILE[$idx]-}"
	local line="${CS_LINE[$idx]-0}"
	local func="${CS_FUNC[$idx]-main}"
	local ts
	ts="$(cs_now)"
	local base="${file##*/}"
	cs_maybe_print_block "$ev" "$base" "$line" "$cmd" "$ec" || return 0
	printf '[%s %s pid=%s %s@%s:%s] FAILED CMD: %s (exit %d)\n' \
		"$ts" "$CS_ROLE" "$CS_PID" "$func" "$base" "$line" "${cmd//$'\n'/\\n}" "$ec" >&2
}

cs_print_wrapper_chain() {
	local printed=0 ts
	ts="$(cs_now)"
	local j=$(((CS_IDX + CS_RING_MAX - 1) % CS_RING_MAX))
	local last=""
	while ((printed < CS_MAX_WRAPPERS)); do
		local w="${CS_CMD[$j]-}"
		[[ -z "$w" ]] && break
		if cs_is_wrapper "$w"; then
			if [[ "$w" != "$last" ]]; then
				local wf="${CS_FILE[$j]-}" wl="${CS_LINE[$j]-0}" fu="${CS_FUNC[$j]-main}"
				printf '[%s %s pid=%s %s@%s:%s] CALLED FROM: %s\n' \
					"$ts" "$CS_ROLE" "$CS_PID" "$fu" "${wf##*/}" "$wl" "${w//$'\n'/\\n}" >&2
				last="$w"
				((printed++))
			fi
			j=$(((j + CS_RING_MAX - 1) % CS_RING_MAX))
		else
			break
		fi
	done
}

cs_print_caller_trace() {
	local i=0 ts
	ts="$(cs_now)"
	while caller "$i" >/dev/null 2>&1; do
		local l f fl
		read -r l f fl < <(caller "$i")
		: "${f:=main}"
		# Inline the source line (best effort)
		local src=""
		if [[ -r "$fl" && "$l" =~ ^[0-9]+$ ]]; then
			src="$(sed -n "${l}p" -- "$fl" 2>/dev/null | tr -d '\r' | tr -s '[:space:]' ' ')"
		fi
		if [[ -n "$src" ]]; then
			printf '[%s %s pid=%s %s@%s:%s] STACK: %s\n' \
				"$ts" "$CS_ROLE" "$CS_PID" "$f" "${fl##*/}" "$l" "$src" >&2
		else
			printf '[%s %s pid=%s %s@%s:%s] STACK\n' \
				"$ts" "$CS_ROLE" "$CS_PID" "$f" "${fl##*/}" "$l" >&2
		fi
		i=$((i + 1))
	done
}

cs_on_err() {
	local ec="${1:-1}"
	local idx
	idx="$(cs_find_relevant_idx)"
	cs_print_failed_at ERR "$idx" "$ec"
	cs_print_wrapper_chain
	cs_print_caller_trace
}

cs_on_exit() {
	local st="$?"
	((st == 0)) && return 0
	local idx
	idx="$(cs_find_relevant_idx)"
	cs_print_failed_at EXIT "$idx" "$st"
	cs_print_wrapper_chain
	cs_print_caller_trace
}

cs_enable() {
	set -o errtrace
	set -o functrace
	trap 'cs_on_debug' DEBUG
	trap 'ec=$?; trap - DEBUG; cs_on_err "$ec"; trap "cs_on_debug" DEBUG' ERR
	trap 'cs_on_exit' EXIT
}
cs_disable() { trap - DEBUG ERR EXIT; }
