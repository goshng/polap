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
[[ -n "${!_POLAP_INCLUDE_}" ]] && return 0
declare "$_POLAP_INCLUDE_=1"
#
################################################################################

source "$script_dir/polap-constants.sh"

################################################################################
# Function to convert base pairs to the highest appropriate unit
# Example usage
# bp=31846726397
# convert_bp $bp
################################################################################
function _polap_utility_convert_bp() {
	local bp=$1

	if ((bp >= 1000000000)); then
		echo "$(bc <<<"scale=1; $bp/1000000000") Gbp"
	elif ((bp >= 1000000)); then
		echo "$(bc <<<"scale=1; $bp/1000000") Mbp"
	elif ((bp >= 1000)); then
		echo "$(bc <<<"scale=1; $bp/1000") kbp"
	else
		echo "$bp bp"
	fi
}

################################################################################
# Function to check if commands are available, taking an array as argument
################################################################################
function check_commands() {
	local cmd_array=("$@") # Capture all the passed arguments into an array
	for cmd in "${cmd_array[@]}"; do
		command -v "$cmd" >/dev/null 2>&1 || {
			echo >&2 "$cmd: not installed"
			return $RETURN_FAIL
		}
	done
	return $RETURN_SUCCESS
}

###############################################################################
# Checks if required main commands are available.
# called early in the code such as reset menu.
###############################################################################
function run_check1() {
	local commands=(
		"bc"
		"seqkit"
		"minimap2"
		"flye"
		"makeblastdb"
		"tblastn"
		"bedtools"
		"prefetch"
		"jellyfish"
		"csvtk"
	)

	# Pass the array elements to the check_commands function
	return $(check_commands "${commands[@]}")
}

###############################################################################
# Checks if FMLRC related commands are available.
# called by prepare-polishing menu.
###############################################################################
function run_check2() {
	local commands=(
		"msbwt"
		"ropebwt2"
		"fmlrc"
	)

	# Pass the array elements to the check_commands function
	return $(check_commands "${commands[@]}")
}

###############################################################################
# Checks if ncbitools related commands are available.
# called by fetch
###############################################################################
function run_check3() {
	local commands=(
		"prefetch"
		"vdb-validate"
		"fasterq-dump"
	)

	# Pass the array elements to the check_commands function
	return $(check_commands "${commands[@]}")
}

###############################################################################
# Checks if ncbitools related commands are available.
###############################################################################
function run_check_ncbitools() {
	local commands=(
		"makeblastdb"
		"tblastn"
		"prefetch"
	)

	# Pass the array elements to the check_commands function
	return $(check_commands "${commands[@]}")
}

###############################################################################
# Function to prompt for confirmation
###############################################################################
function confirm() {
	while true; do
		read -p "$1 [y/n]: " response
		case "$response" in
		[Yy]*) return 0 ;;
		[Nn]*) return 1 ;;
		*) echoerr "Please answer yes (y) or no (n)." ;;
		esac
	done
}

###############################################################################
# Ouputs error messages for executing polap without a proper conda environment.
###############################################################################
function error_polap_conda() {
	_polap_log0 "ERROR: change your conda environment to polap."
	_polap_log0 "SUGGESTION: do the following steps for setting up the polap conda environment"
	_polap_log0 "$ git clone https://github.com/goshng/polap.git"
	_polap_log0 "$ curl -OL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
	_polap_log0 "$ bash Miniconda3-latest-Linux-x86_64.sh"
	_polap_log0 "(base) $ conda update -y -n base conda"
	_polap_log0 "(base) $ conda config --prepend channels bioconda"
	_polap_log0 "(base) $ conda config --prepend channels conda-forge"
	_polap_log0 "(base) $ conda env create -f src/environment.yaml"
	_polap_log0 "(base) $ conda activate polap"
}

# Helper function for checking if a file exists
function check_file_existence() {
	local file=$1
	if [ ! -s "$file" ]; then
		die "ERROR: No such file: $file"
		exit $EXIT_ERROR
	fi
}

function check_folder_existence() {
	local folder=$1
	if [ ! -d "$folder" ]; then
		die "ERROR: No such folder: $folder"
		exit $EXIT_ERROR
	fi
}
