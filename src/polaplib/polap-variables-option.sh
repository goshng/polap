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
# Polap has many global variables, most of which are defined in
# 1. polap-variables-common.sh
# 2. polap-variables-mtcontigs.sh
# This script helps to see all of the global variables and their assigned
# values.
################################################################################

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
	echo "[ERROR] This script must be sourced, not executed: use 'source $BASH_SOURCE'" >&2
	return 1 2>/dev/null || exit 1
fi
: "${_POLAP_DEBUG:=0}"
: "${_POLAP_RELEASE:=0}"

local _arg_outdir="o"
local _arg_long_reads="l.fq"
local _arg_inum="0"
local _arg_jnum="1"
local _arg_single_min="3000"
local _arg_reference=""
local _arg_type="pt"
local _arg_gfa="assemble_graph.gfa"

# Argument parsing
while [[ $# -gt 0 ]]; do
	case "$1" in
	-o)
		_arg_outdir="$2"
		shift 2
		;;
	-i)
		_arg_inum="$2"
		shift 2
		;;
	-j)
		_arg_jnum="$2"
		shift 2
		;;
	-l)
		_arg_long_reads="$2"
		shift 2
		;;
	-w)
		_arg_single_min="$2"
		shift 2
		;;
	-t)
		_arg_type="$2"
		shift 2
		;;
	--reference)
		_arg_reference="$2"
		shift 2
		;;
	--gfa)
		_arg_gfa="$2"
		shift 2
		;;
	-h | --help)
		_polap_log0 "Usage: $FUNCNAME -o _arg_outdir -l _arg_long_reads --reference GFA -w 3000"
		return 0
		;;
	*)
		_polap_log0 "[ERROR] $FUNCNAME Unknown option: $1"
		return 1
		;;
	esac
done
