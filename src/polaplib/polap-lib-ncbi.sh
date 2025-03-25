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

# 1. Using esearch and esummary to Query Genome Size
# You can fetch genome size information from NCBI's Assembly database:
_polap_lib_ncbi-query-genome-size() {
	local _species="${1}"

	echo $(esearch -db assembly -query "${_species}[Organism]" |
		esummary |
		grep 'total_length' |
		awk -F'[<>]' '{print $3}' |
		awk '{sum+=$1; count+=1} END {if (count > 0) print sum/count}' |
		awk '{print int($1)}')

	# xtract -pattern DocumentSummary -element species,assembly_name,total_length
}
