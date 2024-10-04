source "$script_dir/polap-constants.sh"

# Function to convert base pairs to the highest appropriate unit
# Example usage
# bp=31846726397
# convert_bp $bp
_polap_utility_convert_bp() {
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

# Function to check if commands are available, taking an array as argument
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
# Original function that defines the array and calls check_commands
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

function run_check_ncbitools() {
	local commands=(
		"makeblastdb"
		"tblastn"
		"prefetch"
	)

	# Pass the array elements to the check_commands function
	return $(check_commands "${commands[@]}")
}

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
