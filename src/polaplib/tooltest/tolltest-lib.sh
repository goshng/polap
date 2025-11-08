#!/usr/bin/env bash
# polaplib/tooltest/tooltest-lib.sh
# Minimal test DSL for CLI tools (assertions, normalizers, runners)
# Works standalone; if your polap libs exist, it will use them.
set -Eeuo pipefail

: "${TOOLTEST_OUTDIR:=o/tooltest}" # where artifacts go
: "${TOOLTEST_JUNIT:=}"            # path to JUnit report (optional)
: "${TOOLTEST_VERBOSITY:=1}"       # 0..2 (2 = verbose)

mkdir -p -- "${TOOLTEST_OUTDIR}"

# --- tiny helpers --------------------------------------------------------------
_ts() { date '+%Y-%m-%d %H:%M:%S'; }
_here() { printf '%s:%s' "${BASH_SOURCE[1]##*/}" "${BASH_LINENO[0]}"; }
_envt() { [[ -n "${CONDA_DEFAULT_ENV:-}" ]] && printf ' (%s)' "$CONDA_DEFAULT_ENV"; }
logi() { printf '[%s%s %s] %s\n' "$(_ts)" "$(_envt)" "$(_here)" "$*" >&2; }
logv() { [[ "${TOOLTEST_VERBOSITY}" -ge 2 ]] && logi "$@"; }

# Try to enable your failsafe if it exists (nice crash lines).
if [[ -r "$(dirname "${BASH_SOURCE[0]}")/../polap-lib-failsafe.sh" ]]; then
	# shellcheck disable=SC1091
	source "$(dirname "${BASH_SOURCE[0]}")/../polap-lib-failsafe.sh"
	polap_enable_failsafe || true
fi

# --- assertions ---------------------------------------------------------------
assert_file() { [[ -s "$1" ]] || {
	logi "FAIL: file missing/empty: $1"
	return 1
}; }
assert_nonempty() { [[ -s "$1" ]] || {
	logi "FAIL: empty: $1"
	return 1
}; }
assert_exit_ok() {
	local rc="$1"
	shift
	[[ "$rc" -eq 0 ]] || {
		logi "FAIL: exit=$rc: $*"
		return "$rc"
	}
}

# naive semver compare (3 components)
_semver_cmp() {
	awk -v a="$1" -v b="$2" '
  BEGIN{split(a,x,".");split(b,y,".");for(i=1;i<=3;i++){x[i]+=0;y[i]+=0;
       if(x[i]<y[i]){print -1;exit} if(x[i]>y[i]){print 1;exit}} print 0}'
}
assert_version_range() { # cmd min max
	local patt="$1" min="$2" max="$3"
	local v="$($patt 2>&1 | head -n1 | grep -Eo '[0-9]+(\.[0-9]+)+' || true)"
	[[ -n "$v" ]] || {
		logi "FAIL: could not detect version from: $patt"
		return 1
	}
	logv "version detected: $v (need [$min,$max])"
	(($(_semver_cmp "$v" "$min") >= 0)) || {
		logi "FAIL: $v < $min"
		return 1
	}
	(($(_semver_cmp "$v" "$max") <= 0)) || {
		logi "FAIL: $v > $max"
		return 1
	}
}

# --- normalizers (fight non-determinism) --------------------------------------
norm_sort() { sort -S1% -T "${TOOLTEST_OUTDIR}" -o "$1" "$1"; }
norm_strip_ts() { sed -E -i 's/(date|time)=\S+//g;s/^#.*$//g' "$1"; }
norm_float_dp() { # file dp
	local f="$1" dp="${2:-3}"
	awk -v OFS='\t' -v dp="$dp" '
    function round(x,d){p=10^d;return int(x*p+0.5)/p}
    {for(i=1;i<=NF;i++){if($i ~ /^[0-9]*\.[0-9]+$/){$i=round($i,dp)}}; print}' "$f" >"$f.tmp" && mv "$f.tmp" "$f"
}

# --- golden (normalized) compare ----------------------------------------------
golden_compare() { # produced expected [normalizers...]
	local got="$1" exp="$2"
	shift 2
	for n in "$@"; do "$n" "$got"; done
	if ! diff -u --strip-trailing-cr "$exp" "$got" >&2; then
		logi "FAIL: golden mismatch vs $exp"
		return 1
	fi
}

# --- command runners (standalone) ---------------------------------------------
# Minimal argv runner with stdout to file + a repro script
run_argv_to() { # out -- argv...
	local out="$1"
	shift
	[[ "$1" == "--" ]] && shift
	local dest="${TOOLTEST_OUTDIR}/${out}"
	mkdir -p -- "$(dirname "$dest")"
	# repro script
	local repro="${dest%.out}.cmd.sh"
	repro="${repro%.txt}.cmd.sh"
	printf '#!/usr/bin/env bash\nset -Eeuo pipefail\n' >"$repro"
	printf 'cd "%s"\n' "$PWD" >>"$repro"
	printf '%q' "$@" >>"$repro"
	printf '\n' >>"$repro"
	chmod +x "$repro"
	# preview
	local s=""
	for a in "$@"; do printf -v s '%s %q' "$s" "$a"; done
	s="${s# }"
	logi "RUN → ${s}  > ${dest}"
	# run
	"$@" >"$dest"
	local rc=$?
	[[ $rc -eq 0 ]] || logi "non-zero exit: $rc"
	return $rc
}

# --- optional JUnit XML -------------------------------------------------------
_junit_start=0
junit_begin() {
	[[ -n "$TOOLTEST_JUNIT" ]] || return 0
	: >"$TOOLTEST_JUNIT"
	printf '<testsuite>\n' >>"$TOOLTEST_JUNIT"
	_junit_start=1
}
junit_end() { [[ -n "$TOOLTEST_JUNIT" && $_junit_start -eq 1 ]] && printf '</testsuite>\n' >>"$TOOLTEST_JUNIT"; }
_junit_case() {
	local name="$1" status="$2" msg="$3"
	[[ -n "$TOOLTEST_JUNIT" ]] || return 0
	printf '  <testcase name="%s">' "$name" >>"$TOOLTEST_JUNIT"
	if [[ "$status" != "ok" ]]; then
		printf '<failure message="%s"/>' "$(printf '%s' "$msg" | sed 's/[<&>]/ /g')" >>"$TOOLTEST_JUNIT"
	fi
	printf '</testcase>\n' >>"$TOOLTEST_JUNIT"
}

run_case() { # name bash_fn
	local name="$1" fn="$2"
	logi "CASE: $name"
	if "$fn"; then
		_junit_case "$name" ok ""
		return 0
	else
		_junit_case "$name" fail "$name failed"
		return 1
	fi
}
