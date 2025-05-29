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

# Example usage:

# Set some configuration values.
# _polap_set_config_value config.cfg "username" "alice"
# _polap_set_config_value config.cfg "password" "secret123"
# _polap_set_config_value config.cfg "port" "8080"

# Read configuration values.
# echo "Username: $(_polap_get_config_value config.cfg username)"
# echo "Password: $(_polap_get_config_value config.cfg password)"
# echo "Port: $(_polap_get_config_value config.cfg port)"

_polap_create_config() {
	local config_file="$1"
	touch "$config_file"
}

# Function: _polap_get_config_value
# Reads the value associated with a given key from the configuration file.
# Usage: _polap_get_config_value <config> <key>
_polap_get_config_value() {
	local config_file="$1"
	local key="$2"
	# Use grep to find the key and cut to extract the value
	# The regex ensures we match lines that start with the key followed by '='.
	local value
	value=$(grep -m 1 "^${key}=" "$config_file" | cut -d'=' -f2-)
	echo "$value"
}

# Function: _polap_set_config_value
# Writes (or updates) a key-value pair in the configuration file.
# Usage: _polap_set_config_value <key> <value>
_polap_set_config_value() {
	local config_file="$1"
	local key="$2"
	local value="$3"

	# Check if the key exists in the config file.
	if grep -q "^${key}=" "$config_file"; then
		# If the key exists, update the line using sed.
		sed -i.bak "s/^${key}=.*/${key}=${value}/" "$config_file"
	else
		# If the key does not exist, append the key=value pair to the file.
		echo "${key}=${value}" >>"$config_file"
	fi
}
