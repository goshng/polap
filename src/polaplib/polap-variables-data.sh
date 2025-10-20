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

# Global option variables
#
# _brg_outdir
# _brg_sindex
#
# option:
# opt_t_arg

_brg_outdir="${_brg_outdir%/}"

# if [[ -z "${_brg_title}" ]]; then
# 	if [[ "${FUNCNAME[1]}" == run-* ]]; then
# 		local _brg_title="${FUNCNAME[1]#run-}"
# 	fi
# 	_brg_title="${_brg_title%%_*}"
# fi

# _brg_title: may or may not be set
# FUNCNAME[1]: e.g. "run-ptgaul_genus_species"

local _brg_title
if [[ -z "${_brg_title:-}" ]]; then
	_brg_title="${FUNCNAME[1]}"

	# remove suffix like _genus_species
	_brg_title="${_brg_title%%_*}"

	# remove prefix run- if present
	if [[ "${_brg_title}" == run-* ]]; then
		_brg_title="${_brg_title#run-}"
	fi
fi

local _brg_tmpdir="${_brg_outdir}/tmp"
local _brg_target="${_brg_outdir}-${_brg_sindex}"
local _brg_adir="${opt_t_arg:-t1}"
local subcmd1="${_brg_title}"
local _brg_outdir_t="${_brg_outdir}/${opt_t_arg}"
local _brg_outdir_0="${_brg_outdir_t}/0"
local _brg_outdir_i="${_brg_outdir}/${_brg_adir}/${_brg_sindex}"
local _brg_rundir="${_brg_outdir_i}/${_brg_title}"
local _timing_txt="${_brg_outdir_i}/timing-${_brg_title}.txt"
local _stdout_txt="${_brg_outdir_i}/stdout-${_brg_title}.txt"
local _memlog_file="${_brg_outdir_i}/memlog-${_brg_title}.csv"
local _summary_file="${_brg_outdir_i}/summary-${_brg_title}.txt"
local _brg_threads="$(($(grep -c ^processor /proc/cpuinfo)))"

if [[ -v _long["${_brg_target}"] ]]; then
	local long_sra="${_long["${_brg_target}"]}"
else
	echo "Error: ${_brg_target} is not in the CSV."
	return
fi
