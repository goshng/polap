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
# Functions for subcommand template ...
# Describe what they are and what they do.
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

function _run_bolap_install {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - install tools used by bolap

Synopsis:
  bolap ${bolap_cmd} <TOOL>

  bolap ${bolap_cmd} <TOOL>=0.5

Description:
  Install tools; all, apptainer, bandage, bolap, cflye, conda, dflye, fmlrc, fmlrc2, getorganelle, latex, man, mbg, minimal, mitohifi, novoplasty, nvim, oatk, pmat, polap, tippo

Options:
  -l FASTQ
    reads data file

Examples:
  Install tippo:
    bolap ${bolap_cmd} tippo

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_brg_menu[1]} == "help" || "${_brg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	local item="${_brg_menu[1]}"
	local tool="${item%%=*}"
	local version=""
	[[ "$item" == *"="* ]] && version="${item#*=}"

	install-${tool}_genus_species "$version"

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_bolap_setup {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - setup tools

Synopsis:
  bolap ${bolap_cmd} <TOOL>

Description:
  bolap ${bolap_cmd}

Examples:
  Setup polap:
    bolap ${bolap_cmd} polap

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_brg_menu[1]} == "help" || "${_brg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	local item="${_brg_menu[1]}"

	setup_genus_species "${item}"

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_bolap_update {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - setup tools

Synopsis:
  bolap ${bolap_cmd} dev
  bolap ${bolap_cmd} github
  bolap ${bolap_cmd} local

Description:
  bolap ${bolap_cmd} make the content up-to-date.

Examples:
  Setup polap:
    bolap ${bolap_cmd} polap

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_brg_menu[1]} == "help" || "${_brg_help}" == "on" ]]; then
		if [[ ${_brg_menu[1]} == "help" || ${_brg_menu[1]} == "0" ]]; then
			local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		else
			local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		fi
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	# local item="${_brg_menu[1]}"

	"${bolap_cmd}_genus_species" ${_positionals[@]:1}

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_bolap_run {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - execute tools such as organelle genome assembly pipeline.

Synopsis:
  bolap ${bolap_cmd} TOOL

Description:
  bolap ${bolap_cmd} uses a tool to execute something. TOOL includes:
readassemble,
disassemble,
syncassemble,
getorganelle,
ptgaul,
tippo,
oatk,
pmat,


Examples:
  Execute polap readassemble:
    bolap ${bolap_cmd} readassemble Vigna_radiata

  Execute polap disassemble:
    bolap ${bolap_cmd} disassemble Vigna_radiata

  Execute polap syncassemble:
    bolap ${bolap_cmd} syncassemble Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ "${_brg_help}" == "on" ]]; then
		subcmd1="${_brg_menu[1]}"
		_subcmd1_clean="run_${subcmd1//-/_}"
		declare -n ref="help_message_${_subcmd1_clean}"
		local manfile=$(_polap_lib_man-convert_help_message "$ref" "${_brg_menu[0]}-${_brg_menu[1]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	if [[ ${_brg_menu[1]} == "0" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	# local item="${_brg_menu[1]}"
	# local folder=""
	# if [[ "${_brg_menu[2]}" != "0" ]]; then
	# 	folder="${_brg_menu[2]}"
	# fi
	# "${bolap_cmd}_genus_species" "${item}" "${folder}" "${_brg_menu[3]}"
	"${bolap_cmd}_genus_species" ${_positionals[@]:1}

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_bolap_remove {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - 

Synopsis:
  bolap ${bolap_cmd} 

Description:
  bolap ${bolap_cmd} 

Examples:
  Setup polap:
    bolap ${bolap_cmd}

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_brg_menu[1]} == "help" || "${_brg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	local item="${_brg_menu[1]}"
	local folder=""
	if [[ "${_brg_menu[2]}" != "0" ]]; then
		folder="${_brg_menu[2]}"
	fi

	"${bolap_cmd}_genus_species" "${item}" "${folder}" "${_brg_menu[3]}"

	return 0
}

function _run_bolap_download {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - 

Synopsis:
  bolap ${bolap_cmd} 

Description:
  bolap ${bolap_cmd} 

Examples:
  Setup polap:
    bolap ${bolap_cmd}

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_brg_menu[1]} == "help" || "${_brg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	"${bolap_cmd}_genus_species" ${_positionals[@]:1}

	return 0
}

function _run_bolap_config {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${_brg_menu[1]}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - view, add, info

Synopsis:
  bolap ${bolap_cmd} view [column1,column2,...] [query key]

Description:
  bolap ${bolap_cmd} 

Examples:
  List all:
    bolap ${bolap_cmd} view

  List long and short columns of all:
    bolap config view long,short

  List long and short columns of key with Salix:
    bolap config view long,short Salix

  Info:
    bolap config info

  Add a new record to a.csv at the present directory:
    bolap config add Salix_purpurea 1 SRR21824870 NA ONT a.csv
    bolap config add Salix_purpurea 1 SRR21824870 NA PACBIO_SMRT a.csv

  Add a new field to all rows:
    bolap -c a.csv config add-field fieldname value a.csv 

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_brg_menu[1]} == "0" || "${_brg_help}" == "on" ]]; then
		local manfile=$(_bolap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	local bolap_cmd="config-${_brg_menu[1]}"
	subcmd1="${bolap_cmd}"
	echo "subcmd: ${subcmd1}" 2>&1
	if declare -f "${bolap_cmd}_genus_species" >/dev/null 2>&1; then
		"${bolap_cmd}_genus_species" ${_positionals[@]:2}
	else
		_log_echo0 "No such subcommand for config: ${bolap_cmd}"
	fi

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_bolap_benchmark {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - execute organelle genome assembly pipelines

Synopsis:
  bolap ${bolap_cmd} SPECIES-FOLDER [SPECIES-FOLDER2 ...]

Description:
  bolap ${bolap_cmd} uses Oatk, TIPPo, ptGAUL, GetOrganelle, PMAT to assemble
organelle genomes and profiles computing time and memory.
SPECIES-FOLDER is binary species name delimited by an underscore.

Examples:
  Benchamrk of polap and other pipelines:
    bolap ${bolap_cmd} Carex_pseudochinensis

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# _brg_menu defualt: run 0 0 0 ...

	# Display help message
	if [[ "${_brg_help}" == "on" ]]; then
		subcmd1="${_brg_menu[0]}"
		_subcmd1_clean="${subcmd1//-/_}"
		declare -n ref="help_message_${_bolap_type}_${_subcmd1_clean}"
		local manfile=$(_polap_lib_man-convert_help_message "$ref" "${_brg_menu[0]}-${_bolap_type}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	if [[ ${_brg_menu[1]} == "0" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	subcmd1="${bolap_cmd}"
	"${bolap_cmd}_genus_species" ${_positionals[@]:1}

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_bolap_clean {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - 

Synopsis:
  bolap ${bolap_cmd} 

Description:
  bolap ${bolap_cmd} 

Examples:
  Setup polap:
    bolap ${bolap_cmd}

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_brg_menu[1]} == "help" || "${_brg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	local item="${_brg_menu[1]}"
	local folder=""
	if [[ "${_brg_menu[2]}" != "0" ]]; then
		folder="${_brg_menu[2]}"
	fi

	"${bolap_cmd}_genus_species" "${item}" "${folder}" "${_brg_menu[3]}"

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_bolap_get {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - 

Synopsis:
  bolap ${bolap_cmd} 

Description:
  bolap ${bolap_cmd} 

Examples:
  Setup polap:
    bolap ${bolap_cmd}

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_brg_menu[1]} == "help" || "${_brg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	local item="${_brg_menu[1]}"
	local folder=""
	if [[ "${_brg_menu[2]}" != "0" ]]; then
		folder="${_brg_menu[2]}"
	fi

	"${bolap_cmd}_genus_species" "${item}" "${folder}" "${_brg_menu[3]}"

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_bolap_search {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - search tools

Synopsis:
  bolap ${bolap_cmd} <TOOL> [any|end|start]

Description:
  bolap ${bolap_cmd} lists commands.

Examples:
  Search for install:
    bolap ${bolap_cmd} install start

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_brg_menu[1]} == "help" || "${_brg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	local item="${_brg_menu[1]}"
	local where="any"
	if [[ "${_brg_menu[2]}" != "0" ]]; then
		where="${_brg_menu[2]}"
	fi

	list_genus_species "${item}" "${where}"

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_bolap_archive {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - 

Synopsis:
  bolap ${bolap_cmd} 

Description:
  bolap ${bolap_cmd} 

Examples:
  Setup polap:
    bolap ${bolap_cmd}

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_brg_menu[1]} == "help" || "${_brg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	local item="${_brg_menu[1]}"
	local folder=""
	if [[ "${_brg_menu[2]}" != "0" ]]; then
		folder="${_brg_menu[2]}"
	fi

	"${bolap_cmd}_genus_species" "${item}" "${folder}" "${_brg_menu[3]}"

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_bolap_list {
	_run_bolap_search "$@"
}

function _run_bolap_help {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - 

Synopsis:
  bolap ${bolap_cmd} 

Description:
  bolap ${bolap_cmd} 

Examples:
  Setup polap:
    bolap ${bolap_cmd}

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_brg_menu[1]} == "help" || "${_brg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	local item="${_brg_menu[1]}"
	local folder=""
	if [[ "${_brg_menu[2]}" != "0" ]]; then
		folder="${_brg_menu[2]}"
	fi

	if declare -f "${bolap_cmd}_genus_species" >/dev/null 2>&1; then
		"${bolap_cmd}_genus_species" ${_positionals[@]:1}
	fi

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_bolap_use {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - choose polap's data analysis

Synopsis:
  bolap ${bolap_cmd} DATA

Description:
  bolap ${bolap_cmd} allows to choose polap's data analysis; DATA can be aflye, cflye, read.
Each DATA has its own benchmark and man bolap subcommand.

Examples:
  Choose read polap's data analysis.
    bolap ${bolap_cmd} read

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_brg_menu[1]} == "help" || "${_brg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	local item="${_brg_menu[1]}"
	local folder=""
	if [[ "${_brg_menu[2]}" != "0" ]]; then
		folder="${_brg_menu[2]}"
	fi

	subcmd1=use
	_log_echo "${bolap_cmd}_genus_species" "${item}" "${folder}" "${_brg_menu[3]}"
	echo "${item}" >"${HOME}/.bolaprc"

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_bolap_man {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - 

Synopsis:
  bolap ${bolap_cmd} 

Description:
  bolap ${bolap_cmd} 

Examples:
  Setup polap:
    bolap ${bolap_cmd}

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_brg_menu[1]} == "help" || "${_brg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	"${bolap_cmd}_genus_species" ${_positionals[@]:1}

	# store them
	# rest=("${_positionals[@]:1}")
	# echo "${_positionals[@]}"

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_bolap_test {
	local _brg_outdir="${_brg_menu[1]}"
	local _brg_sindex="${_brg_menu[2]:-0}"

	local bolap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - execute tools such as organelle genome assembly pipeline.

Synopsis:
  bolap ${bolap_cmd} TOOL

Description:
  bolap ${bolap_cmd} uses a tool to execute something. TOOL includes:
readassemble,
disassemble,
syncassemble,
getorganelle,
ptgaul,
tippo,
oatk,
pmat,


Examples:
  Execute polap readassemble:
    bolap ${bolap_cmd} readassemble Vigna_radiata

  Execute polap disassemble:
    bolap ${bolap_cmd} disassemble Vigna_radiata

  Execute polap syncassemble:
    bolap ${bolap_cmd} syncassemble Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ "${_brg_help}" == "on" ]]; then
		subcmd1="${_brg_menu[1]}"
		_subcmd1_clean="run_${subcmd1//-/_}"
		declare -n ref="help_message_${_subcmd1_clean}"
		local manfile=$(_polap_lib_man-convert_help_message "$ref" "${_brg_menu[0]}-${_brg_menu[1]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	if [[ ${_brg_menu[1]} == "0" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
		man "$manfile"
		rm -f "$manfile"
		return
	fi

	# local item="${_brg_menu[1]}"
	# local folder=""
	# if [[ "${_brg_menu[2]}" != "0" ]]; then
	# 	folder="${_brg_menu[2]}"
	# fi
	# "${bolap_cmd}_genus_species" "${item}" "${folder}" "${_brg_menu[3]}"
	"${bolap_cmd}_genus_species" ${_positionals[@]:1}

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}

function _run_bolap {

	local bolap_cmd="${_brg_menu[0]}"
	local subcmd1="${_brg_menu[0]}"
	help_message=$(
		cat <<EOF
Name:
  bolap ${bolap_cmd} - execute tools such as organelle genome assembly pipeline.

Synopsis:
  bolap ${bolap_cmd} TOOL

Description:
  bolap ${bolap_cmd} uses a tool to execute something. TOOL includes:
readassemble,
disassemble,
syncassemble,
getorganelle,
ptgaul,
tippo,
oatk,
pmat,


Examples:
  Execute polap readassemble:
    bolap ${bolap_cmd} readassemble Vigna_radiata

  Execute polap disassemble:
    bolap ${bolap_cmd} disassemble Vigna_radiata

  Execute polap syncassemble:
    bolap ${bolap_cmd} syncassemble Vigna_radiata

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ "${_brg_help}" == "on" || ${_brg_menu[1]} == "0" ]]; then
		subcmd1="${_brg_menu[0]}"
		_subcmd1_clean="${subcmd1//-/_}"
		declare -n ref="help_message_${_subcmd1_clean}"
		if [[ -v ref ]]; then
			local manfile=$(_bolap_lib_man-convert_help_message "$ref" "${_brg_menu[0]}")
			man "$manfile"
			rm -f "$manfile"
		else
			_log_echo0 "No such menu: $subcmd1"
		fi

		return
	fi

	if declare -f "${bolap_cmd}_genus_species" >/dev/null 2>&1; then
		"${bolap_cmd}_genus_species" ${_positionals[@]:1}
	fi

	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
