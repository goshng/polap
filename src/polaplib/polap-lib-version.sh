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

source "${_POLAPLIB_DIR}/polap-git-hash-version.sh"
_polap_version=v0.5.1.1-"${_polap_git_hash_version}"

################################################################################
# Function to convert base pairs to the highest appropriate unit
# Example usage
# bp=31846726397
# convert_bp $bp
################################################################################
function _polap_lib_version {
	echo "$(basename $0) ${_polap_version}"
}

_polap_lib_version-check_flye_version() {
	local min_version="2.9.6"
	local version_output
	version_output=$(flye --version 2>/dev/null) || {
		echo "Flye is not installed or not in PATH."
		return 1
	}

	# Extract "2.9.6" from "2.9.6-b1802"
	local current_version
	current_version=$(echo "$version_output" | grep -oE '^[0-9]+\.[0-9]+\.[0-9]+')

	# Compare versions using sort -V
	if [[ "$(printf '%s\n' "$min_version" "$current_version" | sort -V | head -n1)" == "$min_version" ]]; then
		return 0 # current_version >= min_version
	else
		echo "Flye version too old: found $current_version, required ≥ $min_version"
		return 1
	fi
}

# Returns 0 if $1 >= $2 (ge), else 1.
polap_ver_ge() {
	local A="$1" B="$2"
	# Prefer coreutils sort -V if present
	if command -v sort >/dev/null 2>&1; then
		[[ "$(printf '%s\n' "$B" "$A" | sort -V | head -n1)" == "$B" ]]
		return $?
	fi
	# Fallback: awk numeric field compare
	awk -v A="$A" -v B="$B" '
    function splitv(s,a){ n=split(s,a,/[^0-9]+/); for(i=1;i<=n;i++) if(a[i]=="") a[i]=0; return n }
    BEGIN{
      na=splitv(A,aa); nb=splitv(B,bb); n=(na>nb?na:nb)
      for(i=1;i<=n;i++){
        a=(i in aa?aa[i]+0:0); b=(i in bb?bb[i]+0:0)
        if(a>b){exit 0} if(a<b){exit 1}
      }
      exit 0
    }'
}
