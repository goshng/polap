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
# Tip!
# How to extract commands that were expected:
# src/polap.sh reduce-data --redo -v -v -v 2>&1 | grep -E "^rm|^seqkit|^ln"
################################################################################

################################################################################
# Ensure that the current script is sourced only once
source "$script_dir/run-polap-function-include.sh"
_POLAP_INCLUDE_=$(_polap_include "${BASH_SOURCE[0]}")
set +u
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
set -u
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

source "$script_dir/run-polap-function-utilities.sh"

local _polap_version=0.3.7

print_version_history() {
	local _message=$(
		cat <<HEREDOC
POLAP - Plant organelle DNA long-read assembly pipeline.
version ${_polap_version}

----------
v0.3.8
----------
- the 3rd-party software info in the log
- Add proper column names in annotation tables

----------
v0.3.6
----------
- Option --flye-asm-coverage is added so that --coverage is used only for POLAP
- Option --random-seed is added so that random sampling can be done. Previous versions
- Semi-automatic seed contig selection

----------
v0.2.6
----------
- Bioconda package is available.
HEREDOC
	)

	echo "${_message}"
}

function _run_polap_version() { # display version of all
	_log_command_versions
}
