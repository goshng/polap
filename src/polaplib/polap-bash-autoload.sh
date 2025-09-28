#!/usr/bin/env bash
# polap-bash-autoload.sh (with trace)
set -euo pipefail

# Global (exported) list of sourced files for introspection
: "${POLAP_AUTOLOAD_TRACE:=0}"
POLAP_AUTOLOADED_FILES=()

polap_autoload() {
	local base="${1:-}"
	[[ -z "$base" ]] && {
		echo "[autoload] missing POLAPLIB base" >&2
		return 2
	}
	[[ -d "$base" ]] || {
		echo "[autoload] not a directory: $base" >&2
		return 2
	}

	local libdir="$base/lib" cmddir="$base/cmd" has_new_layout=0
	[[ -d "$libdir" || -d "$cmddir" ]] && has_new_layout=1

	_polap_source_file() {
		local f="$1"
		[[ -f "$f" ]] || return 0
		# shellcheck source=/dev/null
		if [[ "$POLAP_AUTOLOAD_TRACE" == "1" ]]; then
			echo sourcing "$f" ... >&2
		fi
		source "$f"
		POLAP_AUTOLOADED_FILES+=("$f")
		if [[ "$POLAP_AUTOLOAD_TRACE" == "1" ]]; then
			echo "[autoload] sourced: $f" >&2
		fi
	}

	local must_first=("polap-lib-version.sh" "polap-constants.sh" "polap-lib-debug.sh" "polap-lib-log.sh")
	for mf in "${must_first[@]}"; do
		if [[ -f "$libdir/$mf" ]]; then
			_polap_source_file "$libdir/$mf"
		elif [[ -f "$base/$mf" ]]; then
			_polap_source_file "$base/$mf"
		else
			echo "[autoload] FATAL: missing required lib: $mf" >&2
			return 3
		fi
	done

	if [[ $has_new_layout -eq 1 && -d "$libdir" ]]; then
		local f
		for f in "$libdir"/*.sh; do
			[[ -e "$f" ]] || continue
			case "$(basename "$f")" in
			polap-lib-version.sh | polap-constants.sh | polap-lib-debug.sh | polap-lib-log.sh) continue ;;
			esac
			_polap_source_file "$f"
		done
	else
		local f
		for f in "$base"/polap-lib-*.sh; do
			[[ -e "$f" ]] || continue
			case "$(basename "$f")" in
			polap-lib-version.sh | polap-constants.sh | polap-lib-debug.sh | polap-lib-log.sh) continue ;;
			esac
			_polap_source_file "$f"
		done
	fi

	if [[ $has_new_layout -eq 1 && -d "$cmddir" ]]; then
		local f
		for f in "$cmddir"/*.sh; do [[ -e "$f" ]] && _polap_source_file "$f"; done
	else
		local f
		for f in "$base"/polap-cmd-*.sh; do [[ -e "$f" ]] && _polap_source_file "$f"; done
	fi
}
export -f polap_autoload
