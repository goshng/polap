#!/usr/bin/env bash
# polap-bash-bolap-autoload.sh (with trace)
set -euo pipefail

: "${POLAP_AUTOLOAD_TRACE:=0}"
BOLAP_AUTOLOADED_FILES=()

bolap_autoload() {
	local base="${1:-}"
	local bolap_type="${2:-read}"
	[[ -z "$base" ]] && {
		echo "[bolap/autoload] missing POLAPLIB base" >&2
		return 2
	}
	[[ -d "$base" ]] || {
		echo "[bolap/autoload] not a directory: $base" >&2
		return 2
	}

	local libdir="$base/lib" cmddir="$base/cmd" has_new_layout=0
	[[ -d "$libdir" || -d "$cmddir" ]] && has_new_layout=1

	_bolap_src() {
		local f="$1"
		[[ -f "$f" ]] || return 0
		source "$f" # shellcheck source=/dev/null
		BOLAP_AUTOLOADED_FILES+=("$f")
		[[ "$POLAP_AUTOLOAD_TRACE" == "1" ]] && echo "[bolap/autoload] sourced: $f" >&2
	}

	local must_first=("polap-lib-version.sh" "bolap-parsing.sh" "polap-lib-man.sh" "polap-lib-tools.sh" "polap-lib-data.sh")
	for mf in "${must_first[@]}"; do
		if [[ -f "$libdir/$mf" ]]; then
			_bolap_src "$libdir/$mf"
		elif [[ -f "$base/$mf" ]]; then
			_bolap_src "$base/$mf"
		else
			echo "[bolap/autoload] FATAL: missing required lib: $mf" >&2
			return 3
		fi
	done

	if [[ $has_new_layout -eq 1 && -d "$libdir" ]]; then
		local f
		for f in "$libdir"/*.sh; do
			[[ -e "$f" ]] || continue
			case "$(basename "$f")" in
			polap-lib-version.sh | bolap-parsing.sh | polap-lib-man.sh | polap-lib-tools.sh | polap-lib-data.sh) continue ;;
			esac
			_bolap_src "$f"
		done
	else
		local f
		for f in "$base"/polap-lib-*.sh; do
			[[ -e "$f" ]] || continue
			case "$(basename "$f")" in
			polap-lib-version.sh | bolap-parsing.sh | polap-lib-man.sh | polap-lib-tools.sh | polap-lib-data.sh) continue ;;
			esac
			echo _bolap_src "$f"
			_bolap_src "$f"
		done
	fi

	if [[ $has_new_layout -eq 1 && -d "$cmddir" ]]; then
		local f
		for f in "$cmddir"/bolap-cmd-*.sh; do [[ -e "$f" ]] && _bolap_src "$f"; done
	else
		local f
		for f in "$base"/bolap-cmd-*.sh; do [[ -e "$f" ]] && _bolap_src "$f"; done
	fi

	local datasrc=""
	if [[ -f "$base/polap-data-${bolap_type}.sh" ]]; then
		datasrc="$base/polap-data-${bolap_type}.sh"
	elif [[ -f "$libdir/polap-data-${bolap_type}.sh" ]]; then
		datasrc="$libdir/polap-data-${bolap_type}.sh"
	else
		echo "[bolap/autoload] ERROR: data file not found: polap-data-${bolap_type}.sh" >&2
		return 4
	fi
	_bolap_src "$datasrc"
}
export -f bolap_autoload
