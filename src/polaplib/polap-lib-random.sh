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

# global variables in polap-lib-random.sh
_polap_var_random_seed=${_arg_random_seed:-$RANDOM}
_polap_var_random_counter=0
_polap_var_random_number=0

# Example usage
# Uncomment the following line to initialize with a specific seed
# init_random 67890  # Set a new seed and restart the sequence
#
# for i in {1..10}; do
#     echo "Random number $i: $(get_random)"
# done

# Counter file to keep track of sequence position
# _polap_var_COUNTER_FILE="${_arg_outdir}/tmp/random_counter.txt"
# _polap_var_SEED_FILE="${_arg_outdir}/tmp/random_seed.txt"

# _polap_var_random_seed=12345
# _polap_var_random_counter=0

# Function to initialize the seed (default: 12345)
_polap_lib_random-init() {
	local seed=${1:-12345} # Use provided seed or default
	_polap_var_random_seed="$seed"
}

# Function to generate repeatable random numbers
_polap_lib_random-get() {
	_polap_var_random_number=$(awk -v seed="$_polap_var_random_seed" -v counter="$_polap_var_random_counter" \
		'BEGIN { srand(seed); for (i=0; i<=counter; i++) num=int(rand()*10000); print num }')

	# echo $((index + 1)) >"$_polap_var_COUNTER_FILE" # Update the counter
	((_polap_var_random_counter++))
}
