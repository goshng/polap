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
v0.5.2
------
- Add readassemble subcommand for ptDNA/mtDNA assembly

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

###############################################################################
# Logs all commands
###############################################################################
function _log_command_versions {
	local commands=(
		# main
		"bc"
		"seqkit"
		"minimap2"
		"flye"
		"makeblastdb"
		"tblastn"
		"bedtools"
		"prefetch"
		"jellyfish"
		"genomescope2"
		"csvtk"
		# fmlrc polishing
		"msbwt"
		"ropebwt2"
		"fmlrc"
		# sratools
		"prefetch"
		"vdb-validate"
		"fasterq-dump"
		# ncbitools
		"makeblastdb"
		"tblastn"
		"prefetch"
	)

	_polap_log1 "------------------------"
	_polap_log1 "conda environment: polap"
	_polap_log1 "------------------------"
	_polap_log1 "version: minimap2: $(minimap2 --version)"
	_polap_log1 "version: flye: $(flye --version)"
	_polap_log1 "version: bedtools: $(bedtools --version)"
	_polap_log1 "version: jellyfish: $(jellyfish --version)"
	_polap_log1 "version: genomescope2: $(genomescope2 --version)"
	_polap_log1 "version: tblastn: $(tblastn -version | head -1)"
	_polap_log1 "version: makeblastdb: $(makeblastdb -version | head -1)"
	_polap_log1 "version: fasterq-dump: $(fasterq-dump --version | tail -2 | head -1)"
	_polap_log1 "version: vdb-validate: $(vdb-validate --version | tail -2 | head -1)"
	_polap_log1 "version: prefetch: $(prefetch --version | tail -2 | head -1)"
	_polap_log1 "version: seqkit: $(seqkit version)"
	_polap_log1 "version: csvtk: $(csvtk version)"
	_polap_log1 "version: bc: $(bc --version | head -1)"
	_polap_log1 "version: gfatools: $(gfatools version | tail -1)"
	_polap_log1 "version: progressiveMauve: $(progressiveMauve --version 2>&1)"

	_polap_log1 "------------------------------"
	_polap_log1 "conda environment: polap-fmlrc"
	_polap_log1 "------------------------------"
	source $HOME/miniconda3/bin/activate polap-fmlrc
	_polap_log1 "version: msbwt: $(msbwt --version 2>&1)"
	_polap_log1 "version: fmlrc: $(fmlrc -v)"
	_polap_log1 "version: ropebwt2: $(ropebwt2 2>&1 | head -2 | tail -1)"
	conda deactivate
}

function _run_polap_version { # display version of all
	print_version_history
	_log_command_versions
}
