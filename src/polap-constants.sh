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
# Ensure that the current script is sourced only once
source "$script_dir/run-polap-function-include.sh"
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

# Constants
EXIT_SUCCESS=0
EXIT_FAIL=1
EXIT_ERROR=2
RETURN_SUCCESS=0
RETURN_FAIL=1

# More constants
_POLAP_ERR_NO_EDGES_GFA=11
_POLAP_ERR_NO_GENOME_SIZE=12
_POLAP_ERR_NO_NK_FQ=13
_POLAP_ERR_NO_SEEDS=14
_POLAP_ERR_NO_DISK_SPACE=15
_POLAP_ERR_MENU_MAP_READS=51
