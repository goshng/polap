#!/usr/bin/env bash
# polaplib/polap-self-check.sh
# Print environment & tool versions; verify basic read/write.

set -Eeuo pipefail

__PLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=polaplib/polap-lib-versions.sh
source "${__PLIB_DIR}/polap-lib-versions.sh"
# shellcheck source=polaplib/polap-lib-toolcheck.sh
source "${__PLIB_DIR}/polap-lib-toolcheck.sh"

tools_default="minimap2:2.26 iqtree2:2.3 samtools:1.17"
tools="${1:-$tools_default}"

echo "=== polap self-check ==="
echo "Time      : $(date -Is)"
echo "Host      : $(hostname)"
echo "CWD       : $(pwd)"
echo "Conda env : ${CONDA_DEFAULT_ENV:-<none>}"
echo "Bash      : ${BASH_VERSION:-<unknown>}"

pyv="$(python3 --version 2>&1 || true)"
echo "Python3   : ${pyv:-<not found>}"
rv="$(Rscript --version 2>&1 | head -n1 || true)"
echo "Rscript   : ${rv:-<not found>}"
echo "PATH      : $PATH"
echo

# Write test
tmpdir="$(mktemp -d -p "${TMPDIR:-/tmp}" polap.selfcheck.XXXXXX)"
if [[ -d "$tmpdir" ]]; then
	echo "Write test: $tmpdir/ok.txt"
	echo "ok" >"$tmpdir/ok.txt"
	[[ -s "$tmpdir/ok.txt" ]] && echo "Write test: OK" || echo "Write test: FAIL"
	rm -rf -- "$tmpdir"
else
	echo "Write test: cannot create temp dir" >&2
fi
echo

# Tools checks
echo "Tools:"
IFS=' ' read -r -a pairs <<<"$tools"
for p in "${pairs[@]}"; do
	name="${p%%:*}"
	need="${p#*:}"
	need="${need:-0}"
	printf '  - %-10s need >= %s ... ' "$name" "$need"
	if command -v "$name" >/dev/null 2>&1; then
		have=""
		case "$name" in
		minimap2) have="$(polap_tool_version_minimap2 || true)" ;;
		iqtree2) have="$(polap_tool_version_iqtree2 || true)" ;;
		samtools) have="$(polap_tool_version_samtools || true)" ;;
		*) have="$("$name" --version 2>&1 | grep -oE '[0-9]+(\.[0-9]+)+' | head -n1)" ;;
		esac
		if [[ -z "$have" ]]; then
			echo "found, version: <unknown>"
		else
			if polap_ver_ge "$have" "$need"; then
				echo "OK (have $have)"
			else
				echo "FAIL (have $have)"
			fi
		fi
		# Show path(s)
		mapfile -t allp < <(command -v -a "$name" 2>/dev/null | awk '!seen[$0]++')
		for ap in "${allp[@]}"; do
			echo "      -> $ap"
		done
	else
		echo "MISSING"
	fi
done

echo
echo "=== done ==="
