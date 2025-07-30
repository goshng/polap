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
# This script has functions that extract information from text files generated
# by polap scripts.
# Much of the input files are handled by polap's test scripts including
# polap-data-v2.sh or polap-data-cflye script.
# Polap is a bash shell script, which creates output files that are used for
# other bash functions of the script.
# TEST-SCC: net-yet
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

source "${_POLAPLIB_DIR}/polap-lib-unit.sh"
source "${_POLAPLIB_DIR}/polap-lib-timing.sh"

_polap_lib_extract-seqkit_stats() {
	local _fq_stats="$1"

	# Check if input file is readable
	if [[ ! -r "$_fq_stats" ]]; then
		# echo "Error: Log file '$_fq_stats' does not exist or is not readable" >&2
		echo "NA"
		return
	fi

	local _l_sra=$(awk 'NR==2 { sub(/\.fastq$/, "", $1); print $1 }' "${_fq_stats}")
	_l_sra=${_l_sra%_1}

	echo "$_l_sra"
}

_polap_lib_extract-long_sra_size_gb() {
	local _v1="$1"

	# Check that input directory is provided
	if [[ -z "$_v1" ]]; then
		echo "Error: No directory provided to extract_long_sra_size_gb" >&2
		return 1
	fi

	local _long_file="${_v1}/long_total_length.txt"
	local _fallback_file="${_v1}/l.fq.txt"
	local _l_sra_size

	# Check for long_total_length.txt or fallback to l.fq.txt
	if [[ -s "$_long_file" ]]; then
		_l_sra_size=$(<"$_long_file")
	elif [[ -s "$_fallback_file" ]]; then
		_l_sra_size=$(<"$_fallback_file")
	else
		echo "Error: Neither '$_long_file' nor '$_fallback_file' is available or non-empty" >&2
		return 1
	fi

	# Check that value is numeric
	if ! [[ "$_l_sra_size" =~ ^[0-9]+$ ]]; then
		echo "Error: SRA size value is not a valid integer: $_l_sra_size" >&2
		return 1
	fi

	# Ensure the conversion function exists
	if ! command -v _polap_utility_convert_bp &>/dev/null; then
		echo "Error: Conversion function '_polap_utility_convert_bp' not found" >&2
		return 1
	fi

	# Attempt the conversion
	local _l_sra_size_gb
	if ! _l_sra_size_gb=$(_polap_utility_convert_bp "${_l_sra_size}"); then
		echo "Error: _polap_utility_convert_bp failed for input $_l_sra_size" >&2
		return 1
	fi

	echo "$_l_sra_size_gb"
}

_polap_lib_extract-seqkit_stats_sum_len() {
	local _fq_stats="$1"

	# Check if input file is readable
	if [[ ! -r "$_fq_stats" ]]; then
		# echo "Error: Log file '$_fq_stats' does not exist or is not readable" >&2
		echo "NA"
		return
	fi

	awk 'NR==2 { print $5 }' "${_fq_stats}"
}

_polap_lib_extract-genome_size() {
	local _brg_outdir="$1"
	local _genome_size="NA"

	if [[ ! -d "$_brg_outdir" ]]; then
		echo "$_genome_size"
		return
	fi

	if [[ -s "${_brg_outdir}/o/short_expected_genome_size.txt" ]]; then
		_genome_size=$(<"${_brg_outdir}/o/short_expected_genome_size.txt")
	elif [[ -s "${_brg_outdir}/short_expected_genome_size.txt" ]]; then
		_genome_size=$(<"${_brg_outdir}/short_expected_genome_size.txt")
	else
		echo "$_genome_size"
		return
	fi

	_genome_size=${_genome_size%%.*}
	echo "$_genome_size"
}

_polap_lib_extract-time_memory_timing_summary_file() {
	local _v1="$1"
	local _timing_file="${_v1}"

	_memory_gb="NA"
	_total_hours="NA"

	if [[ ! -s "$_timing_file" ]]; then
		# echo "Error: timing file not found or empty: $_timing_file" >&2
		# echo "NA"
		_memory_gb="NA"
		_total_hours="NA"
		return
	fi

	_memory_gb=$(grep 'Net increase' "${_v1}" | awk -F'[()]' '{print $2}' | sed 's/ GB//')
	_total_hours=$(grep 'Elapsed time' "${_v1}" | awk '{print $3}')
	_total_hours=$(_polap_lib_timing-convert_to_hours_or_minutes "${_total_hours}")

	return 0
}

_polap_lib_extract-content() {
	local _v1="$1"

	if [[ ! -s "$_v1" ]]; then
		echo "NA"
		return
	fi

	cat "${_v1}"
}

_polap_lib_extract-time_memory_timing_gnu_time_file() {
	local _v1="$1"
	local _timing_file="${_v1}"

	_memory_gb="NA"
	_total_hours="NA"

	if [[ ! -s "$_timing_file" ]]; then
		echo "Error: timing file not found or empty: $_timing_file" >&2
		return 0
	fi

	local _output
	_output=$(_polap_lib_timing-parse-timing "$_timing_file" 2>/dev/null)

	# validate output: must contain two fields
	if [[ "$_output" =~ ^[^[:space:]]+[[:space:]]+[^[:space:]]+ ]]; then
		read -r _memory_gb _total_hours <<<"$_output"
	else
		echo "Warning: malformed output from _polap_lib_timing-parse-timing" >&2
		echo "  Got: $_output" >&2
	fi

	return 0
}

_polap_lib_extract-fasta_seqlen_plain() {
	local fasta="$1"

	if [[ ! -f "$fasta" ]]; then
		echo -e "0"
		return 0
	fi

	awk '
	/^>/ {
		if (seq != "") {
			print header "\t" length(seq)
		}
		header = substr($0, 2)
		seq = ""
		next
	}
	{
		# Remove all whitespace characters before appending
		gsub(/[[:space:]]/, "", $0)
		seq = seq $0
	}
	END {
		if (seq != "") {
			print length(seq)
		}
	}
	' "$fasta" 2>/dev/null || echo -e "0"

	return 0
}

_polap-lib_extract-short_sra_size_gb() {
	local _v1="$1"
	local _s_sra_size

	if [[ -s "${_v1}/short_total_length.txt" ]]; then
		_s_sra_size=$(<"${_v1}/short_total_length.txt")
	else
		local _s_sra_size1
		local _s_sra_size2
		_s_sra_size1=$(<"${_v1}/s1.fq.txt")
		_s_sra_size2=$(<"${_v1}/s2.fq.txt")
		_s_sra_size=$((_s_sra_size1 + _s_sra_size2))
	fi

	local _s_sra_size_gb
	_s_sra_size_gb=$(_polap_lib_unit-convert_bp "${_s_sra_size}")
	echo "$_s_sra_size_gb"
}

_polap_lib_extract-target_read_coverage() {
	local file="$1"

	if [[ ! -f "$file" ]]; then
		echo "NA"
		return
	fi

	local coverage
	# coverage=$(grep -oP 'target coverage:\s*\K\d+' "$file")
	coverage=$(grep -oP 'target coverage:\s*\K[\d.]+(?=x)' "$file" | tail -1)

	if [[ -n "$coverage" ]]; then
		coverage=$(echo "$coverage" | tr -d '[:space:]') # Remove spaces
		printf "%.1f\n" "$coverage"
	else
		echo "Error: Coverage information not found!"
		return 1
	fi
}

_polap_lib_extract-fasta_seqlen() {
	local fasta="$1"

	# If file does not exist, print length 0 and return
	if [[ ! -f "$fasta" ]]; then
		echo "NA"
		return 0
	fi

	# If bioawk is available, use it
	if command -v bioawk >/dev/null 2>&1; then
		bioawk -c fastx 'NR==1 {print length($seq); exit}' "$fasta" 2>/dev/null || echo "0"
		return 0
	fi

	# If seqkit is available and succeeds
	if command -v seqkit >/dev/null 2>&1; then
		if seqkit fx2tab -l -n "$fasta" | cut -f2 | xargs 2>/dev/null; then
			return 0
		fi
	fi

	# Fallback: extract headers and print 0 length
	_polap_lib_extract-fasta_seqlen_plain "$fasta"
	return 0
}

_polap_lib_extract-fasta_numseq() {
	local fasta="$1"

	# Check if file exists
	if [[ ! -f "$fasta" ]]; then
		echo 0
		return 0
	fi

	# Try seqkit if available
	if command -v seqkit >/dev/null 2>&1; then
		seqkit stats "$fasta" 2>/dev/null | awk 'NR > 1 {print $4}' || echo 0
		return 0
	fi

	# Fallback: count '>' lines with grep
	grep -c '^>' "$fasta" 2>/dev/null || echo 0
	return 0
}

_polap_lib_extract-gfa_numseq() {
	local gfa="$1"

	# Check if file exists
	if [[ ! -f "$gfa" ]]; then
		echo 0
		return 0
	fi

	local TMPFILE=$(mktemp)

	gfatools stat -l ${gfa} >${TMPFILE} 2>/dev/null

	local _summary_gfa_number_segments=$(grep "Number of segments:" "${TMPFILE}" | awk '{print $4}')
	local _summary_gfa_total_segment_length=$(grep "Total segment length:" "${TMPFILE}" | awk '{print $4}')
	local _summary_gfa_number_links=$(grep "Number of links:" "${TMPFILE}" | awk '{print $4}')

	# Fallback: count '>' lines with grep
	echo "${_summary_gfa_number_segments}"
	return 0
}

_polap_lib_extract-gfa_seqlen() {
	local gfa="$1"

	# Check if file exists
	if [[ ! -f "$gfa" ]]; then
		echo 0
		return 0
	fi

	local TMPFILE=$(mktemp)

	gfatools stat -l ${gfa} >${TMPFILE} 2>/dev/null

	local _summary_gfa_number_segments=$(grep "Number of segments:" "${TMPFILE}" | awk '{print $4}')
	local _summary_gfa_total_segment_length=$(grep "Total segment length:" "${TMPFILE}" | awk '{print $4}')
	local _summary_gfa_number_links=$(grep "Number of links:" "${TMPFILE}" | awk '{print $4}')

	# cat "${TMPFILE}" >1.txt

	rm -f "${TMPFILE}"

	# Fallback: count '>' lines with grep
	echo "${_summary_gfa_total_segment_length}"
	return 0
}

_polap_lib_extract-system_info() {
	local info_file="$1"
	local -n ref="$2" # name reference to the associative array

	# Define all expected keys
	local keys=(
		hostname kernel_version os_version cpu_model cpu_cores total_memory
		system_uptime load_average filesystem_type mount_source storage_type
		disk_filesystem disk_size disk_used disk_avail disk_use_percent disk_mount
	)

	# Initialize all values to "NA"
	for key in "${keys[@]}"; do
		ref["$key"]="NA"
	done

	# Exit early if file does not exist
	if [[ ! -f "$info_file" ]]; then
		# echo "[ERROR] File not found: $info_file" >&2
		return
	fi

	while IFS= read -r line; do
		case "$line" in
		"Hostname:"*) ref[hostname]="${line#Hostname: }" ;;
		"Kernel Version:"*) ref[kernel_version]="${line#Kernel Version: }" ;;
		"OS Version:"*) ref[os_version]="${line#OS Version: }" ;;
		"CPU Model:"*)
			model="${line#CPU Model: }"
			ref[cpu_model]="${model#Intel(R) Xeon(R) CPU }"
			;;
		"CPU Cores:"*) ref[cpu_cores]="${line#CPU Cores: }" ;;
		"Total Memory:"*) ref[total_memory]="${line#Total Memory: }" ;;
		"System Uptime:"*) ref[system_uptime]="${line#System Uptime: }" ;;
		"Load Average:"*) ref[load_average]="${line#Load Average: }" ;;
		"Filesystem Type"*) ref[filesystem_type]="${line#Filesystem Type (current dir): }" ;;
		"Mount Source:"*) ref[mount_source]="${line#Mount Source: }" ;;
		*"Storage Type"*) ref[storage_type]="${line#*: }" ;;
		"/dev/"*)
			read -r fs size used avail use mounted <<<"$line"
			ref[disk_filesystem]="$fs"
			ref[disk_size]="$size"
			ref[disk_used]="$used"
			ref[disk_avail]="$avail"
			ref[disk_use_percent]="$use"
			ref[disk_mount]="$mounted"
			;;
		esac
	done <"$info_file"
}

_polap_lib_extract-params() {
	local param_file="$1"
	local -n ref="$2" # name reference to the associative array

	# Define all expected keys
	local keys=(
		I P N R A B M D Alpha Memory
	)

	# Initialize all values to "NA"
	for key in "${keys[@]}"; do
		ref["$key"]="NA"
	done

	# Exit early if file does not exist
	if [[ ! -f "$param_file" ]]; then
		# echo "[ERROR] File not found: $info_file" >&2
		return
	fi

	# Initialize default values
	ref=(
		[I]=-1
		[P]=-1
		[N]=-1
		[R]=-1
		[A]=-1
		[B]=-1
		[M]=-1
		[D]=-1
		[Alpha]=-1
		[Memory]=-1
	)

	# Parse the parameter file
	while IFS=": " read -r key value; do
		case "$key" in
		I | P | N | R | A | B | M | D | Alpha | Memory)
			ref["$key"]="$value"
			;;
		esac
	done <"$param_file"
}

# subsample disassemble result
# parse summary1-ordered.txt
_polap_lib_extract-summary1_ordered() {
	local _summary1_ordered_txt="$1"
	local -n ref="$2" # name reference to the associative array

	# Define all expected keys
	local keys=(
		n1 mode1 sd1 index1 second_line_index
	)

	# Initialize all values to "NA"
	for key in "${keys[@]}"; do
		ref["$key"]="NA"
	done

	# Exit early if file does not exist
	if [[ ! -f "$_summary1_ordered_txt" ]]; then
		# echo "[ERROR] File not found: $info_file" >&2
		return
	fi

	# Extract mode value
	# Extract SD value
	# Extract the first index value
	ref[n1]=$(grep "^#n:" "$_summary1_ordered_txt" | awk 'NR==1 {print $2}')
	ref[mode1]=$(grep "^#mode:" "$_summary1_ordered_txt" | awk '{print $2}')
	ref[sd1]=$(grep "^#sd:" "$_summary1_ordered_txt" | awk '{print $2}')
	ref[index1]=$(grep "^#index:" "$_summary1_ordered_txt" | awk 'NR==1 {print $2}')

	# NOTE: so which one is selected for the next stage?
	# problem: the first sorted and the selected one are different.
	# use either one of the two.
	# Then, we should not check the following unnecessary step of checking.
	# Use #index the first one and do not sort them.
	# Or, sort them and use the first one. Problem is we are not sure which index
	# is used.
	local _output=$(awk -F'\t' 'NR==2 {print $1}' "${_summary1_ordered_txt}")
	read -r _second_line_index <<<"$_output"
	ref[second_line_index]="${_second_line_index}"

	if [[ "${ref[index1]}" != "${_second_line_index}" ]]; then
		echo "ERROR: #index: ${ref[index1]} vs. #2nd line index: ${_second_line_index}" >&2
		echo "See ${_summary1_ordered_txt}" >&2
		echo "------------------------------" >&2
		cat "${_summary1_ordered_txt}" >&2
		echo "------------------------------> skip it!" >&2
		# exit 1
	fi
}

# subsample disassemble result
# parse summary2-ordered.txt
_polap_lib_extract-summary2_ordered() {
	local _summary2_ordered_txt="$1"
	local -n ref="$2" # name reference to the associative array

	# Define all expected keys
	local keys=(
		n2 index2 size rate alpha
	)

	# Initialize all values to "NA"
	for key in "${keys[@]}"; do
		ref["$key"]="NA"
	done

	# Exit early if file does not exist
	if [[ ! -f "$_summary2_ordered_txt" ]]; then
		# echo "[ERROR] File not found: $info_file" >&2
		return
	fi

	ref[n2]=$(grep "^#n:" "$_summary2_ordered_txt" | awk 'NR==1 {print $2}')
	local _output=$(awk -F'\t' 'NR==2 {print $1, $2, $4, $11}' "${_summary2_ordered_txt}")
	local _index2
	local _summary2_rate
	local _summary2_size
	local _summary2_alpha
	read -r _index2 _summary2_size _summary2_rate _summary2_alpha <<<"$_output"
	local _summary2_rate_decimal=$(printf "%.10f" "$_summary2_rate")
	ref[rate]=$(echo "scale=4; $_summary2_rate_decimal / 1" | bc)
	ref[size]=$(_polap_lib_unit-convert_bp "${_summary2_size}")
	ref[alpha]=$(echo "scale=2; $_summary2_alpha / 1" | bc | awk '{printf "%.2f\n", $1}')
}
