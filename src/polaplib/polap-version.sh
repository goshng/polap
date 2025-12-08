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
# Print the version information of polap and its tools.
#
# Example:
# polap version
# TODO: rename: polap-cmd-version.sh
################################################################################

################################################################################
# Tip!
# How to extract commands that were executed:
# src/polap.sh reduce-data --redo -v -v -v 2>&1 | grep -E "^rm|^seqkit|^ln"
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

source "${_POLAPLIB_DIR}/run-polap-function-utilities.sh"

function print_version_history {
	local _message=$(
		cat <<HEREDOC
POLAP - Plant organelle DNA long-read assembly pipeline.
version ${_polap_version}

------
v0.5.3
------
- Add a feature to bolap: Tiara, HiMT-filter + Oatk's syncasm

------
v0.5.1
------
- Add bolap

------
v0.4.3
------
- Add disassemble subcommand for ptDNA assembly

------
v0.3.8
------
- Target: Revision 1

------
v0.3.7
------
- Add proper column names in annotation tables

------
v0.3.6
------
- Add a semi-automatic seed contig selection

------
v0.2.6
------
- Bioconda package is available.

POLAP - Plant organelle DNA long-read assembly pipeline.
version ${_polap_version}
printing out packages employed by POLAP to the log: ${LOG_FILE} ...
HEREDOC
	)

	_polap_log0 "${_message}"
}

function _run_polap_version { # display version of all
	print_version_history
	_log_command_versions
}
