#!/usr/bin/env bash
# polaplib/polap-lib-failsafe.sh
# Version: v0.3.0+env+expand
#
# PURPOSE
#   Reliable crash lines for Bash with exact failing site:
#     [ts (ENV) func@file:LINE] <source line> ; <expanded form>
#   + caller chain (oldest → newest) with each frame’s source line and expansion.
#
# HOW IT WORKS
#   • DEBUG trap snapshots the *next* command’s site (caller 0) into a tiny ring.
#   • ERR trap prints the last snapshot (not the library’s frames), then prints
#     the callers’ source lines.
#
# USAGE (in any script that you want covered)
#   set -Eeuo pipefail
#   set -o errtrace
#   source "${_POLAPLIB_DIR}/polap-lib-failsafe.sh"
#   polap_enable_failsafe
#
# TUNABLES
: "${POLAP_FAILSAFE_MAX_CALLERS:=12}"
: "${POLAP_FAILSAFE_RING:=8}"

# ==== Internal state (tiny ring) ==============================================
# The ring tracks the most recent *real* command + site.
#   FS_CMD[k], FS_FUNC[k], FS_FILE[k], FS_LINE[k]
declare -ag FS_CMD=() FS_FUNC=() FS_FILE=() FS_LINE=()
FS_IDX=0

# basename of this library (to skip our own frames)
_FS_SELF_BASE="$(basename -- "${BASH_SOURCE[0]-$0}")"

# ==== Helpers =================================================================

# conda/mamba env tag e.g. "(polap-polish) "
__fs_env_tag() {
	local e=""
	if [[ -n "${CONDA_DEFAULT_ENV-}" ]]; then
		e="$CONDA_DEFAULT_ENV"
	elif [[ -n "${CONDA_PREFIX-}" ]]; then
		e="$(basename -- "$CONDA_PREFIX")"
	fi
	[[ -n "$e" ]] && printf '(%s) ' "$e" || printf ''
}

# one-line printer in the requested format, with optional expanded form
# args: func file line text [expanded]
_polap__fs_print_one() {
	local func="$1" file="$2" line="$3" text="$4" expanded="${5-}"
	local envt
	envt="$(__fs_env_tag)"
	if [[ -n "$expanded" && "$expanded" != "$text" ]]; then
		printf '[%(%Y-%m-%d %H:%M:%S)T %s%s@%s:%s] %s ; %s\n' -1 \
			"$envt" "$func" "$file" "$line" "$text" "$expanded" >&2
	else
		printf '[%(%Y-%m-%d %H:%M:%S)T %s%s@%s:%s] %s\n' -1 \
			"$envt" "$func" "$file" "$line" "$text" >&2
	fi
}

# trim+collapse whitespace (for showing a single source line)
_polap__fs_oneline() { sed -e 's/\r$//' -e 's/[[:space:]]\+/ /g'; }

# read a source line safely (returns nothing on failure)
# args: file line
_polap__fs_src_of() {
	local f="$1" l="$2"
	[[ -r "$f" && "$l" =~ ^[0-9]+$ && "$l" -ge 1 ]] || return 0
	sed -n "${l}p" -- "$f" 2>/dev/null | _polap__fs_oneline
}

# best-effort expansion of a shell line (DANGEROUS on untrusted input)
# It will perform variable/quote/command-subst via eval echo.
__fs_expand_line() {
	local s="$1"
	local out
	out="$(eval "echo $s" 2>/dev/null || true)"
	printf '%s' "$out"
}

# Should we ignore this frame as "library/noise"?
# args: func filebase
_polap__fs_skip_frame() {
	local f="$1" fb="$2"
	# Skip frames from this library or internal failsafe functions
	[[ "$fb" == "$_FS_SELF_BASE" ]] && return 0
	[[ "$f" == _polap__fs_* ]] && return 0
	return 1
}

# ==== DEBUG snapshot ==========================================================
# Capture the NEXT simple command’s caller 0 site into the ring. Keep it tiny.
_polap__fs_on_debug() {
	# If caller 0 is not available, do nothing.
	local l f fl
	if ! read -r l f fl < <(caller 0); then
		return 0
	fi
	: "${f:=main}"
	local fb="${fl##*/}"

	# Skip failsafe’s own frames
	if _polap__fs_skip_frame "$f" "$fb"; then
		return 0
	fi

	# Record the upcoming command text (best effort; may be empty)
	local cmd="${BASH_COMMAND-}"

	# Save snapshot at FS_IDX
	FS_CMD[$FS_IDX]="$cmd"
	FS_FUNC[$FS_IDX]="$f"
	FS_FILE[$FS_IDX]="$fl"
	FS_LINE[$FS_IDX]="$l"
	FS_IDX=$(((FS_IDX + 1) % POLAP_FAILSAFE_RING))
}

# ==== ERR print ===============================================================
_polap__fs_on_err() {
	# Look at the most recent snapshot (one step behind FS_IDX)
	local head=$(((FS_IDX + POLAP_FAILSAFE_RING - 1) % POLAP_FAILSAFE_RING))
	local f="${FS_FUNC[$head]-}"
	local fl="${FS_FILE[$head]-}"
	local l="${FS_LINE[$head]-}"
	# If we have nothing yet (e.g., early failure), try caller 0 as fallback.
	if [[ -z "$fl" || -z "$l" ]]; then
		local _l _f _fl
		if read -r _l _f _fl < <(caller 0); then
			: "${_f:=main}"
			f="$_f"
			fl="$_fl"
			l="$_l"
		else
			fl="${BASH_SOURCE[1]-${BASH_SOURCE[0]-$0}}"
			l="${BASH_LINENO[0]-${LINENO-0}}"
			f="${FUNCNAME[1]-main}"
		fi
	fi

	# 1) Failing site: print the *source line* at file:line (fallback: BASH_COMMAND),
	#    plus an expanded form for copy/paste.
	local src exp
	src="$(_polap__fs_src_of "$fl" "$l")"
	[[ -z "$src" ]] && src="${BASH_COMMAND-}"
	exp="$(__fs_expand_line "$src")"
	_polap__fs_print_one "$f" "${fl##*/}" "$l" "$src" "$exp"

	# 2) Callers (oldest → newest), skipping our library frames
	local frames=() i=0 count=0
	while ((count < POLAP_FAILSAFE_MAX_CALLERS)) && caller "$i" >/dev/null 2>&1; do
		local L F FL
		read -r L F FL < <(caller "$i")
		: "${F:=main}"
		local FB="${FL##*/}"
		if ! _polap__fs_skip_frame "$F" "$FB"; then
			frames+=("$F" "$FB" "$L")
			((count++))
		fi
		((i++))
	done
	# print in reverse to get oldest → newest after the header
	local n=${#frames[@]}
	if ((n > 0)); then
		local idx=$((n - 3))
		while ((idx >= 0)); do
			local F="${frames[idx]}"
			local FB="${frames[idx + 1]}"
			local L="${frames[idx + 2]}"
			local S E
			S="$(_polap__fs_src_of "$FB" "$L")"
			E="$(__fs_expand_line "$S")"
			_polap__fs_print_one "$F" "$FB" "$L" "${S:-CALLED FROM}" "$E"
			((idx -= 3))
		done
	fi
}

# ==== PUBLIC API ==============================================================
polap_enable_failsafe() {
	set -o errtrace
	set -o functrace
	# DEBUG snapshot for the next command
	trap '_polap__fs_on_debug' DEBUG
	# ERR prints the last snapshot + callers
	trap '_polap__fs_on_err' ERR
}

polap_disable_failsafe() {
	trap - DEBUG ERR
}
