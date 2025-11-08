#!/usr/bin/env bash
# polaplib/tooltest/_lib/common.sh
# Minimal test helpers (no external deps). Bash-safe names (underscores).

set -Eeuo pipefail

# ── include guard ─────────────────────────────────────────────────────────────
if [[ -n "${_POLAP_RUN_COMMON_SOURCED:-}" ]]; then
	return 0
fi

_Polap_Run_Common_Version="v0.1.0"
readonly _Polap_Run_Common_Version
_Polap_Run_Common_File="${BASH_SOURCE[0]:-polap-lib-run-common.sh}"
readonly _Polap_Run_Common_File
_POLAP_RUN_COMMON_SOURCED=1

# ── config knobs ──────────────────────────────────────────────────────────────
: "${POLAP_VERBOSE:=1}" # 0..4 → which levels go to screen (0=quiet)
: "${LOG_FILE:=}"       # if set, all messages also append to this file
: "${POLAP_TAGDATE:=1}" # show date in logs (1=yes)

# Optional: enable your failsafe crash tracer if present
if [[ "${POLAP_FAILSAFE:-0}" == "1" ]]; then
	__TT_THIS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
	__TT_POLAPLIB="$(cd "${__TT_THIS_DIR}/../.." && pwd)"
	if [[ -r "${__TT_POLAPLIB}/polap-lib-failsafe.sh" ]]; then
		# shellcheck disable=SC1091
		source "${__TT_POLAPLIB}/polap-lib-failsafe.sh"
		polap_enable_failsafe
	fi
fi

: "${POLAP_TOOLTEST_OUTDIR:=o/tooltest}"

tt_conda_tag() {
	local tag=""
	if [[ -n "${CONDA_DEFAULT_ENV:-}" ]]; then tag="(${CONDA_DEFAULT_ENV}) "; fi
	printf '%s' "$tag"
}

tt_ts() { date '+%Y-%m-%d %H:%M:%S'; }

tt_banner() {
	printf '[%s %s%s] %s\n' "$(tt_ts)" "$(tt_conda_tag)" "${FUNCNAME[1]:-main}" "$*"
}

tt_esc() {
	local s=""
	local a
	for a in "$@"; do printf -v s '%s %q' "$s" "$a"; done
	printf '%s' "${s# }"
}

tt_echo_cmd() {
	local esc
	esc="$(tt_esc "$@")"
	printf '[%s %s%s] %s ; %s\n' "$(tt_ts)" "$(tt_conda_tag)" "${FUNCNAME[1]:-main}" "$esc" "$*"
}

tt_require_tool() {
	local t="$1"
	if ! command -v "$t" >/dev/null 2>&1; then
		printf '[%s %s%s] SKIP: %s not found in PATH\n' "$(tt_ts)" "$(tt_conda_tag)" "${FUNCNAME[1]:-main}" "$t" >&2
		return 99
	fi
}

tt_assert_rc() {
	local have="$1" want="${2:-0}" msg="${3:-}"
	if [[ "$have" -ne "$want" ]]; then
		printf '[%s %s%s] ASSERT FAIL rc=%d want=%d %s\n' "$(tt_ts)" "$(tt_conda_tag)" "${FUNCNAME[1]:-main}" "$have" "$want" "$msg" >&2
		return 1
	fi
}

tt_assert_file() {
	local f="$1"
	if [[ ! -s "$f" ]]; then
		printf '[%s %s%s] ASSERT FAIL missing/empty: %s\n' "$(tt_ts)" "$(tt_conda_tag)" "${FUNCNAME[1]:-main}" "$f" >&2
		return 1
	fi
}

tt_assert_contains() {
	local pattern="$1" file="$2"
	if ! grep -Fq -- "$pattern" "$file"; then
		printf '[%s %s%s] ASSERT FAIL: %s not found in %s\n' "$(tt_ts)" "$(tt_conda_tag)" "${FUNCNAME[1]:-main}" "$pattern" "$file" >&2
		return 1
	fi
}

tt_mkout() { mkdir -p -- "$1"; }

tt_run() {
	tt_echo_cmd "$@"
	"$@"
}

tt_run_to() {
	local out="$1" err="$2"
	shift 2
	tt_echo_cmd "$@"
	"$@" >"$out" 2>"$err"
}

# Fixtures
tt_ensure_mini_fa() {
	local dst="$1"
	[[ -s "$dst" ]] && return 0
	mkdir -p -- "$(dirname "$dst")"
	cat >"$dst" <<'FA'
>chr1
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC
>chr2
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAA
FA
}

tt_ensure_mini_sam() {
	local dst="$1"
	[[ -s "$dst" ]] && return 0
	mkdir -p -- "$(dirname "$dst")"
	cat >"$dst" <<'SAM'
@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:44
@SQ	SN:chr2	LN:44
r001	0	chr1	1	60	10M	*	0	0	ACGTACGTAA	*	NM:i:0
SAM
}

tt_ensure_mini_aln_fa() {
	local dst="$1"
	[[ -s "$dst" ]] && return 0
	mkdir -p -- "$(dirname "$dst")"
	cat >"$dst" <<'ALN'
>sp1
ACGTACGTACGTACGTACGT
>sp2
ACGTACGTACGTACGTACGA
>sp3
ACGTACGTACGTACGTACGG
>sp4
ACGTACGTACGTACGTACGC
ALN
}

tt_ensure_mini_fq() {
	local dst="$1"
	[[ -s "$dst" ]] && return 0
	mkdir -p -- "$(dirname "$dst")"
	cat >"$dst" <<'FQ'
@r1
ACGTACGTACGTACGTACGTACGTACGTACGTACGT
+
#####################################
@r2
ACGTACGTACGTACGTACGTACGTACGTACGTACGA
+
######################################
FQ
}

# ── small helpers ─────────────────────────────────────────────────────────────
polap__ts() { date '+%Y-%m-%d %H:%M:%S'; }
polap__envtag() {
	local e="${CONDA_DEFAULT_ENV-}"
	[[ -n "$e" ]] && printf '(%s) ' "$e"
}
# tag for the *caller* at stack depth (default = 2: your wrapper calls this)
polap__tag() {
	local depth="${1:-2}"
	local func="${FUNCNAME[$depth]-main}"
	local file="$(basename -- "${BASH_SOURCE[$depth]-${BASH_SOURCE[0]}}")"
	local line="${BASH_LINENO[$((depth - 1))]-0}"
	printf '%s@%s:%s' "$func" "$file" "$line"
}

# join/quote argv for preview (no eval)
polap__preview_argv() {
	local s="" a
	for a in "$@"; do printf -v s '%s %q' "$s" "$a"; done
	printf '%s' "${s# }"
}

# "expand" a string without executing it: only parameter expansion
# (Note: this does *not* run commands; it just expands $vars within the current shell)
polap__expand_string() {
	local s="$*"
	# shellcheck disable=SC2016
	eval "printf '%s' \"$s\""
}

# core emitter: L0..L3 → writes to LOG_FILE (if set) and maybe to screen
# depth: which stack level to print as caller tag (3 works well for wrappers)
polap__emit() {
	local level="$1" depth="$2"
	shift 2
	local msg="$*"
	local ts=""
	[[ "$POLAP_TAGDATE" == "1" ]] && ts="$(polap__ts)"
	local tag="$(polap__tag "$depth")"
	local env="$(polap__envtag)"
	local line
	if [[ -n "$ts" ]]; then
		printf -v line '[%s %s%s] %s' "$ts" "$env" "$tag" "$msg"
	else
		printf -v line '[%s%s] %s' "$env" "$tag" "$msg"
	fi

	# To file
	if [[ -n "$LOG_FILE" ]]; then
		printf '%s\n' "$line" >>"$LOG_FILE"
	fi
	# To screen depending on POLAP_VERBOSE vs level (L0=need 1, L1=2, L2=3, L3=4)
	local need
	case "$level" in
	L0) need=1 ;;
	L1) need=2 ;;
	L2) need=3 ;;
	L3) need=4 ;;
	*) need=2 ;;
	esac
	if ((POLAP_VERBOSE >= need)); then
		printf '%s\n' "$line" >&2
	fi
}

# public log shorthands (depth=3 is suited for wrappers that call emit)
polap_log0() { polap__emit L0 3 "$*"; }
polap_log1() { polap__emit L1 3 "$*"; }
polap_log2() { polap__emit L2 3 "$*"; }
polap_log3() { polap__emit L3 3 "$*"; }

# tiny file helpers
polap__ensure_dir() { install -d -- "$1"; }
polap__mkparent() { install -d -- "$(dirname -- "$1")"; }

# assert helpers (exit 1 on failure)
polap_assert_file() {
	local f="$1"
	[[ -s "$f" ]] || {
		polap__emit L0 3 "ASSERT file missing: $f"
		return 1
	}
}
polap_assert_contains() {
	local needle="$1" hay="$2"
	grep -Fq -- "$needle" "$hay" || {
		polap__emit L0 3 "ASSERT text not found: $needle in $hay"
		return 1
	}
}
