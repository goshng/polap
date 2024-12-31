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

source "$script_dir/polap-constants.sh"
source "$script_dir/run-polap-function-log.sh"

################################################################################
# try a command
################################################################################

function try { "$@" || die "cannot $*"; }

################################################################################
# Function to convert base pairs to the highest appropriate unit
# Example usage
# bp=31846726397
# convert_bp $bp
################################################################################
function _polap_utility_convert_bp {
	local bp=${1%.*}
	if ((bp >= 1000000000000)); then
		echo "$(bc <<<"scale=1; $bp/1000000000000") Tbp"
	elif ((bp >= 1000000000)); then
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
# Function to convert base pairs to the highest appropriate unit
# Example usage
# bp=$(convert_unit_to_bp 1gb)
################################################################################
function _polap_utility_convert_unit_to_bp {
	local input="$1"
	local number unit

	# Separate the numeric part and the unit
	number=$(echo "$input" | grep -oE '^[0-9]+(\.[0-9]+)?')
	unit=$(echo "$input" | grep -oE '[a-zA-Z]+$' | tr '[:upper:]' '[:lower:]')

	# Handle the case where no unit is provided
	if [[ -z "$unit" ]]; then
		echo "$number"
		return 0
	fi

	# Normalize units (handle 'b' and 'bp' equivalently)
	case "$unit" in
	tbp | tb) echo "$(bc <<<"$number * 1000000000000")" ;;
	gbp | gb) echo "$(bc <<<"$number * 1000000000")" ;;
	mbp | mb) echo "$(bc <<<"$number * 1000000")" ;;
	kbp | kb) echo "$(bc <<<"$number * 1000")" ;;
	bp | b) echo "$number" ;;
	*)
		echo "Error: Invalid unit '$unit'" >&2
		return 1
		;;
	esac
}

################################################################################
# Get contig total length
################################################################################
function _polap_utility_get_contig_length {
	local contig_file="$1"
	local length_file="$2"
	_polap_log3_pipe "seqkit stats -Ta $contig_file |
    csvtk cut -t -f sum_len |
    csvtk del-header \
    >${length_file}"
}

################################################################################
# Function to check if commands are available, taking an array as argument
################################################################################
function check_commands {
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
function run_check1 {
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
function run_check2 {
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
function run_check3 {
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
function run_check_ncbitools {
	local commands=(
		"makeblastdb"
		"tblastn"
		"prefetch"
	)

	# Pass the array elements to the check_commands function
	return $(check_commands "${commands[@]}")
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

###############################################################################
# Function to prompt for confirmation
###############################################################################
function confirm {
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
function error_polap_conda {
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
function check_file_existence {
	local file=$1
	if [ ! -s "$file" ]; then
		die "ERROR: No such file: $file"
		exit $EXIT_ERROR
	fi
}

function check_folder_existence {
	local folder=$1
	if [ ! -d "$folder" ]; then
		die "ERROR: No such folder: $folder"
		exit $EXIT_ERROR
	fi
}

# Function to check and display the size of an SRA accession in GB, MB, or KB
function get_sra_size {
	# Check if the user provided an SRA accession
	if [ -z "$1" ]; then
		echo "Usage: get_sra_size <SRA_ACCESSION>"
		return 1
	fi

	SRA_ACCESSION=$1

	# Get the size in bytes using vdb-dump
	size_in_bytes=$(vdb-dump $SRA_ACCESSION --info | grep size | awk -F':' '{print $2}' | tr -d ' ,' | xargs)

	# Check if the size was retrieved
	if [ -z "$size_in_bytes" ]; then
		echo "Error: Could not retrieve size for $SRA_ACCESSION"
		return 1
	fi

	# Convert size to human-readable format
	if [ "$size_in_bytes" -ge 1073741824 ]; then
		size_in_gb=$(echo "scale=2; $size_in_bytes / 1073741824" | bc)
		echo "Size of $SRA_ACCESSION: $size_in_gb GB"
	elif [ "$size_in_bytes" -ge 1048576 ]; then
		size_in_mb=$(echo "scale=2; $size_in_bytes / 1048576" | bc)
		echo "Size of $SRA_ACCESSION: $size_in_mb MB"
	else
		size_in_kb=$(echo "scale=2; $size_in_bytes / 1024" | bc)
		echo "Size of $SRA_ACCESSION: $size_in_kb KB"
	fi
}

################################################################################
# creats a range file
################################################################################
function _create_range {
	# Extract arguments
	local params="$1"
	local output_file="$2"

	# Extract the start, end, and count values from the input string
	IFS=',' read -r start end count <<<"${params}"

	# Initialize array
	local numbers=()

	# Handle the special case where start equals end and count is 1
	if [[ "$start" == "$end" ]]; then
		for ((i = 0; i < count; i++)); do
			numbers+=("$start")
		done
	else
		# Calculate the step
		# Calculate the step size
		local step=$(((end - start) / (count - 1)))

		# Generate numbers and add to array
		for ((i = 0; i < count; i++)); do
			value=$((start + i * step))
			numbers+=("$value")
		done
	fi

	# Save array to a file
	echo "${numbers[@]}" >"$output_file"
}

function _create_range_float {
	# Extract arguments
	local params="$1"
	local output_file="$2"

	# Split params into start, end, and count
	IFS=',' read -r start end count <<<"$params"

	# Initialize array
	local numbers=()

	# Handle the special case where start equals end and count is 1
	if [[ "$start" == "$end" ]]; then
		for ((i = 0; i < count; i++)); do
			numbers+=("$(printf "%.2f" "$start")")
		done
	else
		# Calculate the step
		local step
		step=$(echo "($end - $start) / ($count - 1)" | bc -l)

		# Generate numbers and add to array
		for i in $(seq 0 $((count - 1))); do
			value=$(echo "$start + $i * $step" | bc -l)
			formatted_value=$(printf "%.2f" "$value")
			numbers+=("$formatted_value")
		done
	fi

	# Save the array to the output file
	# printf "%s\n" "${numbers[@]}" >"$output_file"
	echo "${numbers[@]}" >"$output_file"
}

# _polap_tempfiles=() # Initialize an array to hold temp files

# Function to create temp files and add them to the array
function _polap_create_tempfile {
	local tmp=$(mktemp ${_arg_outdir}/tmp/XXX)
	# _polap_tempfiles+=("$tmp")
	echo "$tmp"
}

# Create multiple temp files
# tempfile1=$(create_tempfile)
# tempfile2=$(create_tempfile)

function _polap_delete_tempfiles {
	# Set up a trap to delete all files in the array
	rm -f "${_arg_outdir}/tmp"/*
}

# function _polap_delete_tempfiles {
# 	# Set up a trap to delete all files in the array
# 	trap 'rm -f "${_polap_tempfiles[@]}"' EXIT
# }
#

function _polap_gunzip-fastq {
	# Example file name (replace this with your file or a loop to handle multiple files)
	# fastq_file="example.fastq.gz"
	local fastq_file=$1
	local new_file="${fastq_file}"

	# Check if the file is gzipped
	if file "$fastq_file" | grep -q "gzip compressed"; then
		gunzip "$fastq_file"
		new_file="${fastq_file%.gz}"
	fi
	echo "${new_file}"
}
