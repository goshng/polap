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

# Bash assert, modeled after C's assert()
# Usage:
#   _polap_assert '[[ -s "$file" ]]' "file must be non-empty"
#   _polap_assert '(( x > 0 ))'
#
# Disable all assertions by exporting NDEBUG=1 (like C):
#   export NDEBUG=1
#
# Exit code: 134 (SIGABRT-like). Uses kill -ABRT if available.

# Internal: print a compact stack trace (caller → main)
_polap__trace() {
	local i=1
	# FUNCNAME[0] is _polap__trace, [1] is _polap_assert; start at caller of assert
	while ((i < ${#FUNCNAME[@]})); do
		local fn="${FUNCNAME[$i]:-MAIN}"
		local src="${BASH_SOURCE[$i]:-?}"
		local ln="${BASH_LINENO[$((i - 1))]:-?}"
		printf '    at %s (%s:%s)\n' "$fn" "$src" "$ln" >&2
		((i++))
	done
}

# Internal: abort with SIGABRT if we can; otherwise exit 134
_polap__abort() {
	# If we're in a subshell, kill $$ is safe; otherwise fall back to exit.
	if kill -l ABRT >/dev/null 2>&1; then
		kill -ABRT $$
		# If signal delivery fails for any reason:
		exit 134
	else
		exit 134
	fi
}

# No-op when NDEBUG is set (like C)
if [[ "${NDEBUG:-0}" -eq 1 ]]; then
	_polap_assert() { :; }
else
	_polap_assert() {
		local expr="${1:?_polap_assert: missing expression}"
		shift || true
		# Evaluate the expression in a context that won't be tripped by `set -e`.
		if eval "$expr"; then
			return 0
		fi

		# Assemble error report
		local file="${BASH_SOURCE[1]:-?}"
		local line="${BASH_LINENO[0]:-?}"
		local func="${FUNCNAME[1]:-MAIN}"
		{
			printf 'Assertion failed: (%s)\n' "$expr"
			printf '  Location: %s:%s in %s()\n' "$file" "$line" "$func"
			if (($#)); then printf '  Message: %s\n' "$*"; fi
			_polap__trace
		} >&2

		_polap__abort
	}
fi

# Convenience helpers (optional):
# Assert the previous command succeeded; use after a command:  cmd || _polap_assert_ok "why it must succeed"
_polap_assert_ok() {
	local rc=$?
	((rc == 0)) && return 0
	_polap_assert 'false' "previous command exited with $rc${1:+ – }$*"
}

# Assert two integers equal: _polap_assert_eq "$got" "$want" "numbers must match"
_polap_assert_eq() {
	local got="${1:?}"
	local want="${2:?}"
	shift 2 || true
	_polap_assert "(( ${got} == ${want} ))" "${*:-expected ${want}, got ${got}}"
}
