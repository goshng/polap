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

_polap_lib_extract-seqkit_stats() {
	local _fq_stats="$1"

	# Check if input file is readable
	if [[ ! -r "$_fq_stats" ]]; then
		echo "Error: Log file '$_fq_stats' does not exist or is not readable" >&2
		return 1
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
	awk 'NR==2 { print $5 }' $1
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
		echo "Error: timing file not found or empty: $_timing_file" >&2
		return 0
	fi

	_memory_gb=$(grep 'Net increase' "${_v1}" | awk -F'[()]' '{print $2}' | sed 's/ GB//')
	_total_hours=$(grep 'Elapsed time' "${_v1}" | awk '{print $3}')
	_total_hours=$(_polap_lib_timing-convert_to_hours_or_minutes "${_total_hours}")

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
