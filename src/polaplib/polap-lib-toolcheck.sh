#!/usr/bin/env bash
# polaplib/polap-lib-toolcheck.sh
# Tool presence and minimum-version checks. Friendly messages.

# Resolve our dir and ensure version comparator exists
__PLIB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if ! type polap_ver_ge >/dev/null 2>&1; then
	# shellcheck source=polaplib/polap-lib-versions.sh
	source "${__PLIB_DIR}/polap-lib-version.sh"
fi

# === Tool-specific version extraction (add as needed) =========================

polap_tool_version_minimap2() {
	# minimap2 --version prints e.g. "2.26-r1175"
	local v
	v="$(minimap2 --version 2>/dev/null | head -n1 || true)"
	[[ -n "$v" ]] || return 1
	printf '%s\n' "${v%%-*}"
}

polap_tool_version_iqtree2() {
	iqtree2 -version 2>&1 |
		sed -n 's/^IQ-TREE version \([0-9.]\+\).*/\1/p' |
		head -n1
}

polap_tool_version_samtools() {
	samtools --version 2>&1 |
		sed -n 's/^samtools \([0-9.]\+\).*/\1/p' |
		head -n1
}

polap_tool_version_seqkit() {
	seqkit version 2>&1 |
		sed -n 's/^seqkit v\([0-9.]\+\).*/\1/p' |
		head -n1
}

# === Generic "require tool" ===================================================

polap_require_tool() {
	# usage: polap_require_tool NAME MINVER [HINT]
	local name="$1" need="$2" hint="${3:-install via conda/mamba}"

	if ! command -v "$name" >/dev/null 2>&1; then
		printf '[ERR] missing tool: %s. %s\n' "$name" "$hint" >&2
		exit 3
	fi

	local have=""
	case "$name" in
	minimap2) have="$(polap_tool_version_minimap2 || true)" ;;
	iqtree2) have="$(polap_tool_version_iqtree2 || true)" ;;
	samtools) have="$(polap_tool_version_samtools || true)" ;;
	seqkit) have="$(polap_tool_version_seqkit || true)" ;;
	*)
		# best-effort generic getter
		have="$("$name" --version 2>&1 | grep -oE '[0-9]+(\.[0-9]+)+' | head -n1)"
		;;
	esac

	if [[ -z "$have" ]]; then
		printf '[ERR] cannot parse %s version; need >= %s\n' "$name" "$need" >&2
		exit 4
	fi

	if ! polap_ver_ge "$have" "$need"; then
		printf '[ERR] %s %s found, need >= %s\n' "$name" "$have" "$need" >&2
		exit 4
	fi

	# Optional: record into o/versions.txt if LOG_FILE is defined
	if [[ -n "${LOG_FILE:-}" ]]; then
		printf '%s\t%s\n' "$name" "$have" >>"${LOG_FILE%/*}/versions.txt"
	fi
}
