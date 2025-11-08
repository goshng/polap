#!/usr/bin/env bash
# polaplib/tooltest/run-one.sh
set -Eeuo pipefail
tool="${1:?tool name (dir under tooltest/)}"
tdir="$(cd "$(dirname "$0")" && pwd)/${tool}"
[[ -d "$tdir" ]] || {
	echo "no such tool dir: $tdir" >&2
	exit 2
}

st=0
# Execute any of {smoke,io,gold,perf}.sh if present and executable
for t in smoke io gold perf; do
	f="${tdir}/${t}.sh"
	[[ -x "$f" ]] || continue
	bash "$f" || st=$((st | 1))
done
exit "$st"
