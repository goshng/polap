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
# Convert help_message to unix man page.
#
# manfile=$(convert_help_to_man "$help_message")
# man "$manfile"
# rm -f "$manfile"
################################################################################
_polap_help_message=$(
	cat <<'EOF'
Name:
  gzip, gunzip, zcat - compress or expand files

Synopsis:
  gzip [ -acdfhklLnNrtvV19 ] [-S suffix] [ name ... ]
  gunzip [ -acfhklLnNrtvV ] [-S suffix] [ name ... ]
  zcat [ -fhLV ] [ name ... ]

Description:
  gzip reduces the size of files using Lempel-Ziv (LZ77) compression.
  Files are replaced with a .gz version. gunzip and zcat decompress files.

Options:
  -a, --ascii
    Convert text end-of-lines for non-Unix systems.

  -c, --stdout, --to-stdout
    Write output to stdout; keep original files.

  -d, --decompress, --uncompress
    Decompress instead of compress.

  -f, --force
    Force overwrite of existing files or links.

  -l, --list
    Display compression stats:
    With --verbose:

  -N, --name
    Store or restore original file name and timestamp.

  -S .suf, --suffix .suf
    Use suffix .suf instead of .gz

Advanced usage:
  Concatenate:
    gzip -c file1 > foo.gz
    gzip -c file2 >> foo.gz
    gunzip -c foo.gz

  Recompress:
    gzip -cd old.gz | gzip > new.gz

Environment:
  GZIP environment variable sets default options (deprecated)

Diagnostics:
  Exit status:
    0  success
    1  error
    2  warning

  Examples:
    file: not in gzip format
    file already exists; do you wish to overwrite (y or n)?
    gunzip: corrupt input

Bugs:
  Files >4 GB report incorrect size via --list
  In rare cases, --best compresses worse than default

Copyright:
  Copyright © 1992–1993 Jean-loup Gailly
  Free Software Foundation (1998–2018)

Author:
  Jean-loup Gailly and Mark Adler
EOF
)

function _polap_lib_man-convert_help_message() {
	local help="$1"
	local cmd="polap $2"
	local version="$_polap_version"
	local date=$(date +"%B %Y")
	local tmpfile
	tmpfile=$(mktemp "/tmp/${cmd}.XXXX.1")

	{
		echo ".TH ${cmd^^} 1 \"$date\" \"$cmd $version\" \"General Commands Manual\""

		local section=""
		local in_option_block=0

		while IFS= read -r raw_line; do
			# If line starts with exactly 2 spaces, strip them; otherwise, leave as-is
			if [[ "$raw_line" =~ ^[[:space:]]{2}[^[:space:]] ]]; then
				line="${raw_line:2}"
			else
				line="$raw_line"
			fi

			case "$line" in
			"Name:"*)
				section="NAME"
				echo ".SH NAME"
				continue
				;;
			"Synopsis:"*)
				section="SYNOPSIS"
				echo ".SH SYNOPSIS"
				continue
				;;
			"Description:"*)
				section="DESCRIPTION"
				echo ".SH DESCRIPTION"
				continue
				;;
			"Options:"*)
				section="OPTIONS"
				echo ".SH OPTIONS"
				continue
				;;
			"Advanced usage:"*)
				section="ADVANCED USAGE"
				echo ".SH ADVANCED USAGE"
				continue
				;;
			"Environment:"*)
				section="ENVIRONMENT"
				echo ".SH ENVIRONMENT"
				continue
				;;
			"Diagnostics:"*)
				section="DIAGNOSTICS"
				echo ".SH DIAGNOSTICS"
				continue
				;;
			"Examples:"*)
				section="EXAMPLES"
				echo ".SH EXAMPLES"
				continue
				;;
			"Inputs:"*)
				section="INPUTS"
				echo ".SH INPUTS"
				continue
				;;
			"Outputs:"*)
				section="OUTPUTS"
				echo ".SH OUTPUTS"
				continue
				;;
			"Menus:"*)
				section="MENUS"
				echo ".SH MENUS"
				continue
				;;
			"View:"*)
				section="VIEW"
				echo ".SH VIEW"
				continue
				;;
			"TODO:"*)
				section="TODO"
				echo ".SH TODO"
				continue
				;;
			"Bugs:"*)
				section="BUGS"
				echo ".SH BUGS"
				continue
				;;
			"Copyright:"*)
				section="COPYRIGHT"
				echo ".SH COPYRIGHT"
				continue
				;;
			"See also:"*)
				section="SEE ALSO"
				echo ".SH SEE ALSO"
				continue
				;;
			"Author:"*)
				section="AUTHOR"
				echo ".SH AUTHOR"
				continue
				;;
			"")
				in_option_block=0
				echo
				continue
				;;
			esac

			case "$section" in
			NAME | DESCRIPTION | ADVANCED\ USAGE | ENVIRONMENT | DIAGNOSTICS | EXAMPLES | INPUTS | OUTPUTS | MENUS | VIEW | TODO | BUGS | COPYRIGHT | SEE\ ALSO | AUTHOR)
				echo "$line"
				;;
			SYNOPSIS)
				echo ".B $line"
				echo ".br"
				;;
			OPTIONS)
				if [[ "$line" =~ ^- ]]; then
					echo ".TP"
					echo ".B $line"
					in_option_block=1
				elif [[ "$in_option_block" -eq 1 && "$line" =~ ^[[:space:]]{2,} ]]; then
					echo "       $line"
				elif [[ "$in_option_block" -eq 1 ]]; then
					echo "$line"
				fi
				;;
			esac
		done <<<"$help"

	} >"$tmpfile"

	echo "$tmpfile"
}
