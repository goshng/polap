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

################################################################################
# _polap_lib_lines-skip1 file.txt       # 2..EOF (skip header)
# _polap_lib_lines-skip1 file.txt 1     # whole file (include header)
# _polap_lib_lines-skip1 file.txt 10    # 2..10 (skip header)
# _polap_lib_lines-skip1 file.txt 5 10  # 5..10
# _polap_lib_lines-skip1 file.txt -5    # 2..(last-5)
# _polap_lib_lines-skip1 file.txt -10 -5 # 10th-last..5th-last
# _polap_lib_lines-skip1 file.txt 3 -1  # 3..(last-1)
#
# cat file.txt | _polap_lib_lines-skip1     # 2..EOF
#
# cat file.txt | _polap_lib_lines-skip1 1   # whole file
#
_polap_lib_lines-skip1() {
	local file A B
	A=""
	B=""

	# Parse args
	if [[ $# -ge 1 && (-e "$1" || "$1" == "-") ]]; then
		case $# in
		1) file="$1" ;;
		2)
			file="$1"
			B="$2"
			;;
		3)
			file="$1"
			A="$2"
			B="$3"
			;;
		*)
			echo "Usage: select_lines <file|-> [start] [end]" >&2
			return 1
			;;
		esac
	else
		file="-"
		case $# in
		0) ;;
		1) B="$1" ;;
		2)
			A="$1"
			B="$2"
			;;
		*)
			echo "Usage: select_lines <file|-> [start] [end]" >&2
			return 1
			;;
		esac
	fi

	# Stdin buffering (enables negative indices)
	local tmp=
	if [[ "$file" == "-" ]]; then
		tmp="$(mktemp)"
		cat >"$tmp"
		file="$tmp"
	fi
	_cleanup() { [[ -n "$tmp" && -f "$tmp" ]] && rm -f "$tmp"; }

	if [[ ! -f "$file" ]]; then
		echo "ERROR: file not found: $file" >&2
		_cleanup
		return 1
	fi

	local total_lines
	total_lines=$(wc -l <"$file")
	((total_lines == 0)) && {
		_cleanup
		return
	}

	# Defaults
	if [[ -z "$A" && -z "$B" ]]; then
		A=2
	elif [[ -z "$A" && -n "$B" ]]; then
		if [[ "$B" == "1" ]]; then
			A=1
			B=$total_lines # whole file
		else
			A=2 # up to B, skipping header
		fi
	fi

	_is_int() { [[ "$1" =~ ^-?[0-9]+$ ]]; }

	# Normalize A
	if [[ -n "$A" ]]; then
		if ! _is_int "$A"; then
			echo "ERROR: start must be integer (got: $A)" >&2
			_cleanup
			return 2
		fi
		if [[ "$A" == -* ]]; then
			A=$((total_lines + A + 1))
		fi
		((A < 1)) && A=1
		((A > total_lines)) && {
			_cleanup
			return
		}
	fi

	# Normalize B
	if [[ -n "$B" ]]; then
		if ! _is_int "$B"; then
			echo "ERROR: end must be integer (got: $B)" >&2
			_cleanup
			return 2
		fi
		if [[ "$B" == -* ]]; then
			B=$((total_lines + B + 1))
		fi
		((B < 1)) && {
			_cleanup
			return
		}
		((B > total_lines)) && B=$total_lines

		# NEW: strict check for B < A
		if ((A > B)); then
			echo "ERROR: end (B=$B) is before start (A=$A)" >&2
			_cleanup
			return 3
		fi
	fi

	# Output
	if [[ -n "$B" ]]; then
		sed -n "${A},${B}p" -- "$file"
	else
		tail -n +"$A" -- "$file"
	fi

	_cleanup
}
