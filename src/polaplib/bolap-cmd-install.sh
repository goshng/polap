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
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_brg_menu[0]}")
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

	"${bolap_cmd}_genus_species" "${item}" "${folder}" "${_brg_menu[3]}"

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
