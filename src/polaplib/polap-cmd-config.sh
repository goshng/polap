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

function _run_polap_config {
	# Enable debugging if _POLAP_DEBUG is set
	[ "$_POLAP_DEBUG" -eq 1 ] && set -x
	_polap_log_function "Function start: $(echo $FUNCNAME | sed s/_run_polap_//)"

	# Set verbosity level: stderr if verbose >= 2, otherwise discard output
	local _polap_output_dest="/dev/null"
	[ "${_arg_verbose}" -ge "${_polap_var_function_verbose}" ] && _polap_output_dest="/dev/stderr"

	# Grouped file path declarations
	source "${_POLAPLIB_DIR}/polap-variables-common.sh" # '.' means 'source'

	local polap_cmd="${FUNCNAME##*_}"
	help_message=$(
		cat <<EOF
Name:
  polap ${polap_cmd} - create, load, and update polap config

Synopsis:
  polap ${polap_cmd} [options]

Description:
  A set of config is saved as a YAML config file in $HOME/.polap/profiles.

Options:
  -l FASTQ
    reads data file

Examples:
  Create default profiles:
    polap ${polap_cmd} init

  List profiles:
    polap ${polap_cmd} list

  Copy profile source to dest:
    polap ${polap_cmd} copy source dest
    polap ${polap_cmd} cp source dest

  View profile hifi:
    polap ${polap_cmd} view hifi

  Edit profile hifi:
    polap ${polap_cmd} [edit] hifi --preset hifi [options]

Copyright:
  Copyright © 2025 Sang Chul Choi
  Free Software Foundation (1998–2018)

Author:
  Sang Chul Choi
EOF
	)

	# Display help message
	if [[ ${_arg_menu[1]} == "help" || "${_arg_help}" == "on" ]]; then
		local manfile=$(_polap_lib_man-convert_help_message "$help_message" "${_arg_menu[0]}")
		man "$manfile" >&3
		rm -f "$manfile"
		return
	fi

	# Display the content of output files
	if [[ "${_arg_menu[1]}" == "view" ]]; then
		local preset="base"

		if [[ "${_arg_menu[2]}" != "outfile" ]]; then
			preset="${_arg_menu[2]}"
		fi

		python3 "${_POLAPLIB_DIR}/polap-py-config-validate.py" \
			--preset "${preset}" \
			--print >&3

		_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
		# Disable debugging if previously enabled
		[ "$_POLAP_DEBUG" -eq 1 ] && set +x
		return 0
		exit $EXIT_SUCCESS
	fi

	if [[ "${_arg_menu[1]}" == init ]]; then
		python3 "${_POLAPLIB_DIR}/polap-py-config-template.py" --force
		cp -p hifi.yaml "${HOME}/.polap/profiles"
		cp -p ont.yaml "${HOME}/.polap/profiles"
	fi

	if [[ "${_arg_menu[1]}" == "new" ]]; then
		source "${_POLAPLIB_DIR}/polap-bash-config-new.sh"

		polap_config_new \
			--preset mt_hifi \
			--template hifi \
			--config-dir "${HOME}/.polap/profiles" \
			--wgs-mode \
			--reads /data/all.fq.gz \
			--anchors /data/anchors/mt.id.all.txt \
			--hmm-db /data/hmm/mt_genes.hmm \
			--outdir /work/autotune \
			--threads 16 \
			--k 121 --s 27 --no-hpc \
			--final-if "genes_score < 0.95 || breadth < 0.97" \
			--c-radius "0,10,20" \
			--min-shared 5 --jaccard-min 0.015 --topk-nei 50 --steps 2
	fi

	if [[ "${_arg_menu[1]}" == "load" ]]; then
		local PCFG_PRESET="mt_hifi"
		source "${_POLAPLIB_DIR}/polap-bash-config-load.sh"

		_polap_log0 "$PCFG_PRESET"   # hifi
		_polap_log0 "$PCFG_WGS_MODE" # true/false
		_polap_log0 "$PCFG_READS"    # /data/all.fq.gz
	fi

	if [[ "${_arg_menu[1]}" == "list" ]]; then
		_polap_log0 "Configs: ${HOME}/.polap/profiles"
		ls -1 "${HOME}/.polap/profiles" >&3
	fi

	if [[ "${_arg_menu[1]}" == "copy" || "${_arg_menu[1]}" == "cp" ]]; then
		_polap_log0 "copy config preset:${_arg_menu[2]} to preset:${_arg_menu[3]}"

		cp -p "${HOME}/.polap/profiles/${_arg_menu[2]}.yaml" \
			"${HOME}/.polap/profiles/${_arg_menu[3]}.yaml"
	fi

	if [[ "${_arg_menu[1]}" == "infile" || "${_arg_menu[1]}" == "edit" ]]; then
		_polap_log0 "update the config:${_arg_config_path}"
	fi

	_polap_log3 "Function end: $(echo $FUNCNAME | sed s/_run_polap_//)"
	# Disable debugging if previously enabled
	[ "$_POLAP_DEBUG" -eq 1 ] && set +x
	return 0
}
