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
# Config the polap outdir.
# TODO: rename: polap-lib-config.sh
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

# polap_pcfg_dump_yaml: dump PCFG_* env vars to YAML
# Usage:
#   polap_pcfg_dump_yaml --out config.yaml
#   polap_pcfg_dump_yaml --prefix PCFG > config.yaml
polap_pcfg_dump_yaml() {
	local out="" prefix="PCFG"
	while (($#)); do
		case "$1" in
		--out)
			out="${2:?}"
			shift 2
			;;
		--prefix)
			prefix="${2:?}"
			shift 2
			;;
		-h | --help)
			cat <<EOF
Usage: polap_pcfg_dump_yaml [--out FILE] [--prefix NAME]
  --out FILE     write YAML to FILE (default: stdout)
  --prefix NAME  read NAME_* envs (default: PCFG_)

Examples:
  export PCFG_OUTDIR=o; export PCFG_SPECIES='Macadamia tetraphylla'
  polap_pcfg_dump_yaml --out config.yaml
EOF
			return 0
			;;
		*)
			echo "Unknown option: $1" >&2
			return 2
			;;
		esac
	done

	local pfx="${prefix}_"
	# Collect matching env var names deterministically
	mapfile -t __vars < <(compgen -v | awk -v pfx="$pfx" 'index($0,pfx)==1' | sort -f)

	# Writer: stdout or file
	if [[ -n "$out" ]]; then
		: >"$out" || {
			echo "ERROR: cannot write $out" >&2
			return 1
		}
	fi
	__emit() {
		if [[ -n "$out" ]]; then printf '%s\n' "$*" >>"$out"; else printf '%s\n' "$*"; fi
	}

	# Iterate and render YAML
	local v key key_lc val low
	for v in "${__vars[@]}"; do
		# Key: strip prefix and lowercase
		key="${v#${pfx}}"
		key_lc="$(tr '[:upper:]' '[:lower:]' <<<"$key")"

		# Read value safely even under set -u
		if [[ -z "${!v+x}" ]]; then
			val=""
		else
			val="${!v}"
		fi

		# Trim leading/trailing spaces
		val="${val#"${val%%[![:space:]]*}"}"
		val="${val%"${val##*[![:space:]]}"}"
		low="${val,,}"

		# Multiline -> block scalar
		if [[ "$val" == *$'\n'* ]]; then
			__emit "${key_lc}: |"
			while IFS= read -r line; do __emit "  ${line}"; done <<<"$val"
			continue
		fi

		# List support: comma-separated -> sequence
		if [[ "$val" == *,* ]]; then
			__emit "${key_lc}:"
			IFS=',' read -r -a __arr <<<"$val"
			for item in "${__arr[@]}"; do
				# trim
				item="${item#"${item%%[![:space:]]*}"}"
				item="${item%"${item##*[![:space:]]}"}"
				ilow="${item,,}"
				if [[ "$item" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then
					__emit "  - $item"
				elif [[ "$ilow" =~ ^(on|true|yes|1)$ ]]; then
					__emit "  - true"
				elif [[ "$ilow" =~ ^(off|false|no|0)$ ]]; then
					__emit "  - false"
				else
					item="${item//\'/\'\'}"
					__emit "  - '${item}'"
				fi
			done
			continue
		fi

		# Scalar rendering
		if [[ "$val" =~ ^-?[0-9]+([.][0-9]+)?$ ]]; then
			__emit "${key_lc}: ${val}"
		elif [[ "$low" =~ ^(on|true|yes|1)$ ]]; then
			__emit "${key_lc}: true"
		elif [[ "$low" =~ ^(off|false|no|0)$ ]]; then
			__emit "${key_lc}: false"
		else
			# quote strings safely
			val="${val//\'/\'\'}"
			__emit "${key_lc}: '${val}'"
		fi
	done
}
