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

# Find the best *.organelle.chloroplast.gfa in a given directory.
# Selection: highest roundN; tie-breaker: highest filter.N
# Prints the chosen file path on stdout, or "NONE" if not found.
#
# chosen=$(find_best_chloro_gfa /tmp/emptydir) || true
# echo "chosen=[$chosen]"
#
_polap_lib_tools-find_best_chloro_gfa() {
	local dir="${1:?Usage: find_best_chloro_gfa <directory>}"
	[[ -d "$dir" ]] || {
		echo "ERROR: not a directory: $dir" >&2
		echo "NONE"
		return 2
	}

	local best="" best_round=-1 best_filter=-1
	local matched=false

	for f in "$dir"/*organelle.chloroplast.gfa; do
		if [[ -e "$f" ]]; then
			matched=true
		else
			continue
		fi

		local base r=0 filt=-1
		base=$(basename -- "$f")

		# roundN (missing -> 0)
		if [[ $base =~ (^|[.])round([0-9]+)([.]|$) ]]; then
			r=${BASH_REMATCH[2]}
		fi
		# filter.N (missing -> -1)
		if [[ $base =~ (^|[.])filter\.([0-9]+)([.]|$) ]]; then
			filt=${BASH_REMATCH[2]}
		fi

		if ((r > best_round || (r == best_round && filt > best_filter))); then
			best="$f"
			best_round=$r
			best_filter=$filt
		fi
	done

	if ! $matched; then
		echo "NONE"
		return 1
	fi

	echo "$best"
}
