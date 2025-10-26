#!/usr/bin/env bash
# polap-bash-bolap-autoload-min.sh
# Only load the allow-listed libs for bolap, then bolap commands, then data sheet.

set -euo pipefail
: "${POLAP_AUTOLOAD_TRACE:=0}"
BOLAP_AUTOLOADED_FILES=()

bolap_autoload_min() {
	local base="${1:-}" bolap_type="${2:-read}"
	[[ -z "$base" ]] && {
		echo "[bolap/autoload-min] missing POLAPLIB base" >&2
		return 2
	}
	[[ -d "$base" ]] || {
		echo "[bolap/autoload-min] not a directory: $base" >&2
		return 2
	}

	local libdir="$base/lib" cmddir="$base/cmd" new_layout=0
	[[ -d "$libdir" || -d "$cmddir" ]] && new_layout=1

	_src() {
		local f="$1"
		[[ -f "$f" ]] || return 0
		# shellcheck source=/dev/null
		source "$f"
		BOLAP_AUTOLOADED_FILES+=("$f")
		if [[ "$POLAP_AUTOLOAD_TRACE" == "1" ]]; then
			echo "[bolap/autoload] sourced: $f" >&2
		fi
	}

	# 1) Allow-list of libraries (exactly what you listed)
	local ALLOW_LIBS=(
		"polap-lib-version.sh"
		"polap-lib-conda.sh"
		"polap-lib-tools.sh"
		"polap-lib-timing.sh"
		"polap-lib-unit.sh"
		"polap-lib-array.sh"
		"polap-lib-number.sh"
		"polap-lib-data.sh"
		"polap-lib-file.sh"
		"polap-lib-process.sh"
		"polap-lib-extract.sh"
		"polap-lib-csv.sh"
		"polap-lib-man.sh"
	)

	# Helper to resolve base or libdir
	_resolve_and_source() {
		local name="$1"
		if [[ $new_layout -eq 1 && -f "$libdir/$name" ]]; then
			_src "$libdir/$name"
		elif [[ -f "$base/$name" ]]; then
			_src "$base/$name"
		else
			echo "[bolap/autoload-min] WARN: missing lib (skipped): $name" >&2
		fi
	}

	local L
	for L in "${ALLOW_LIBS[@]}"; do
		_resolve_and_source "$L"
	done

	# 2) Commands (UI first, then engine commands if you want them available)
	if [[ $new_layout -eq 1 && -d "$cmddir" ]]; then
		local f
		for f in "$cmddir"/bolap-cmd-*.sh; do [[ -e "$f" ]] && _src "$f"; done
	else
		local f
		for f in "$base"/bolap-cmd-*.sh; do [[ -e "$f" ]] && _src "$f"; done
	fi

	# 3) Data sheet
	local datasrc=""
	if [[ -f "$base/polap-data-${bolap_type}.sh" ]]; then
		datasrc="$base/polap-data-${bolap_type}.sh"
	elif [[ -f "$libdir/polap-data-${bolap_type}.sh" ]]; then
		datasrc="$libdir/polap-data-${bolap_type}.sh"
	else
		echo "[bolap/autoload-min] ERROR: data file not found: polap-data-${bolap_type}.sh" >&2
		return 4
	fi
	_src "$datasrc"
}
export -f bolap_autoload_min
