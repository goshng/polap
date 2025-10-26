################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

################################################################################
# Convert numbers between different units.
################################################################################

################################################################################
# Ensure that the current script is sourced only once
source "${_POLAPLIB_DIR}/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u
if [[ -n "${!_POLAP_INCLUDE_}" ]]; then
	set -u
	return 0
fi
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
	echo "[ERROR] This script must be sourced, not executed: use 'source $BASH_SOURCE'" >&2
	return 1 2>/dev/null || exit 1
fi
: "${_POLAP_DEBUG:=0}"
: "${_POLAP_RELEASE:=0}"

# Version: v0.4.0
# _polap_lib_math-percentify VALUE [DECIMALS] [ROUNDING]
# - VALUE: numeric in [0,1] (supports scientific notation, e.g. 1.23e-2). REQUIRED.
# - DECIMALS: non-negative integer, default 2.
# - ROUNDING: nearest|floor|ceil  (default: nearest)
#
# Prints: percentage string (e.g., "87.35") and exits 0.
# On invalid input: prints error with caller info, exits 1.

_polap_lib_math-percentify() {
	local value="$1"
	local decimals="${2:-2}"
	local mode="${3:-nearest}"

	# --- caller info for error reporting ---
	local _src="${BASH_SOURCE[1]:-$(basename "$0")}"
	local _line="${BASH_LINENO[0]:-?}"
	local _func="${FUNCNAME[1]:-}" # may be empty if called from "main"
	local _where="$_src:$_line"
	[[ -n "$_func" && "$_func" != "source" ]] && _where+=" in ${_func}()"

	# --- validate decimals ---
	if ! [[ "$decimals" =~ ^[0-9]+$ ]]; then
		printf '[percentify] Invalid DECIMALS "%s" (need non-negative integer). Called from %s.\n' "$decimals" "$_where" >&2
		exit 1
	fi

	# --- validate rounding mode ---
	case "$mode" in
	nearest | floor | ceil) ;;
	*)
		printf '[percentify] Invalid ROUNDING "%s" (use nearest|floor|ceil). Called from %s.\n' "$mode" "$_where" >&2
		exit 1
		;;
	esac

	# --- validate numeric (no negatives), allow sci-notation ---
	# Forms: 123, 123., 123.45, .45, 1e-3, 1.2E+4, etc. Optional leading '+', but NOT '-'.
	local num_re='^[+]?(([0-9]+([.][0-9]*)?|[.][0-9]+))([eE][+-]?[0-9]+)?$'
	if ! [[ "$value" =~ $num_re ]]; then
		printf '[percentify] Invalid VALUE "%s" (not a number in [0,1]). Called from %s.\n' "$value" "$_where" >&2
		exit 1
	fi

	# --- bound check: must be 0 <= value <= 1 ---
	# Use awk for safe numeric compare (handles sci-notation).
	if ! awk -v x="$value" 'BEGIN{exit !(x>=0 && x<=1)}'; then
		printf '[percentify] VALUE "%s" out of range [0,1]. Called from %s.\n' "$value" "$_where" >&2
		exit 1
	fi

	# --- compute percentage with rounding ---
	# Use awk for math & rounding (handles sci-notation). Non-negative domain simplifies floor/ceil logic.
	case "$mode" in
	nearest)
		awk -v x="$value" -v n="$decimals" 'BEGIN{r=x*100; printf "%.*f", n, r}'
		;;
	floor)
		awk -v x="$value" -v n="$decimals" 'BEGIN{p=1; for(i=0;i<n;i++) p*=10; r=x*100; y=int(r*p)/p; printf "%.*f", n, y}'
		;;
	ceil)
		awk -v x="$value" -v n="$decimals" 'BEGIN{p=1; for(i=0;i<n;i++) p*=10; r=x*100; t=r*p; y=(t==int(t)?t:int(t)+1)/p; printf "%.*f", n, y}'
		;;
	esac
}
