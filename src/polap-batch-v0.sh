#!/bin/bash

_polap_script_bin_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" || {
	echo "Couldn't determine the script's running directory, which probably matters, bailing out" >&2
	exit 2
}

_polap_data_cmd="${_polap_script_bin_dir}/polap-data-v0.sh"

echo "A: " ${_polap_data_cmd} test test
${_polap_data_cmd} test test
