#!/usr/bin/env bash
# polaplib/polap-lib-run.sh
# Helpers: preview commands (quoted & unquoted), safe run, allowed codes,
# atomic write to final path, and simple assertions.

polap_preview_expanded() {
	# Print both a quoted (paste-safe) and an unquoted (quick-typed) form.
	# usage: polap_preview_expanded cmd arg1 arg2 ...
	local quoted="" a
	for a in "$@"; do printf -v quoted '%s %q' "$quoted" "$a"; done
	local unquoted
	unquoted=$(printf '%s ' "$@")
	printf 'quoted  :%s\n' "${quoted# }"
	printf 'unquoted:%s\n' "${unquoted% }"
}

polap_run() {
	# usage: polap_run cmd arg...
	local ts envtag
	ts="$(date '+%Y-%m-%d %H:%M:%S')"
	envtag="${CONDA_DEFAULT_ENV:+ (${CONDA_DEFAULT_ENV})}"
	printf '[%s%s RUN] ' "$ts" "${envtag:-}" >&2
	polap_preview_expanded "$@" >&2
	"$@"
}

polap_run_ok() {
	# usage: polap_run_ok "0 1" cmd arg...  (treat listed exit codes as success)
	local ok="$1"
	shift
	"$@"
	local st=$?
	for c in $ok; do [[ $st -eq $c ]] && return 0; done
	return $st
}

polap_atomic_write() {
	# usage: polap_atomic_write FINAL_PATH -- cmd ...
	local final="$1"
	shift
	[[ "$1" == "--" ]] && shift
	local tmp dir
	dir="$(dirname -- "$final")"
	install -d -- "$dir"
	tmp="$(mktemp -p "${TMPDIR:-/tmp}" ".polap.tmp.XXXXXX")"
	if "$@" >"$tmp"; then
		install -D -m 0644 -- "$tmp" "$final"
		rm -f -- "$tmp"
	else
		local st=$?
		echo "[ERR] producer failed for $final (exit $st)" >&2
		rm -f -- "$tmp"
		return $st
	fi
}

polap_assert_file() {
	# usage: polap_assert_file /path/to/file
	[[ -r "$1" ]] || {
		echo "[ERR] no such readable file: $1" >&2
		exit 5
	}
}

polap_assert_dir_writable() {
	# usage: polap_assert_dir_writable /path/to/dir
	install -d -m 0755 -- "$1" 2>/dev/null || {
		echo "[ERR] directory not writable: $1" >&2
		exit 6
	}
}
