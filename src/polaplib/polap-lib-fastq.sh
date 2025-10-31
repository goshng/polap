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
# Bash functions for handling FASTQ files.
# TEST-SCC: not yet
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

# arg1: input fasta or fastq file
# arg2: output text file with the total number of bases in the input sequence
# file
_polap_lib_fastq-total-length-of() {
	local input1="$1"
	local output="$2"

	check_file_existence "${input1}"

	_polap_log3_pipe "seqkit stats -Ta ${input1} -j 4 |\
		csvtk cut -t -f sum_len |\
		csvtk del-header \
		>${output}"

	return 0
}

# Function to handle FASTQ files
# output: fastq file
# input1: fastq file
# input2: fastq file (optional)
#
# Example usage
# _polap_concatenate_fastq_to "s1.fq" "file1.fq"
# _polap_concatenate_fastq_to "s1.fq" "file1.fq" "file2.fq.gz"
# _polap_concatenate_fastq_to "s1.fq" "file1.fq.gz" "file2.fq.gz"
_polap_lib_fastq-concatenate-fastq-files() {
	local output_file="$1"     # Output file (e.g., s1.fq)
	local input_file1="$2"     # First input file
	local input_file2="${3:-}" # Second input file (optional)

	# Check if the first input file exists
	if [[ ! -f "$input_file1" ]]; then
		die "ERROR: no such file: $input_file1"
	fi

	# Create or overwrite the output file
	rm -f "${output_file}"

	# If only one file is provided and it is not compressed, create a soft link
	if [[ -z "$input_file2" && ! "$input_file1" == *.gz ]]; then
		ln -s "$(realpath "$input_file1")" "$output_file" || {
			die "ERROR: cannot create a soft link."
		}
		_polap_log2 "  soft link created: $output_file -> $input_file1"
		return 0
	fi

	# Function to process a single file
	(
		process_file() {
			local file="$1"
			if [[ $file == *.gz ]]; then
				gunzip -c "$file" >>"$output_file" || {
					die "ERROR: gunzip -c $file"
				}
			else
				cat "$file" >>"$output_file" || {
					die "ERROR: cat $file"
				}
			fi
		}

		# Process the first file
		process_file "$input_file1"

		# Process the second file if provided
		if [[ -n "$input_file2" ]]; then
			if [[ ! -f "$input_file2" ]]; then
				die "ERROR: no such file: $input_file2"
			fi
			process_file "$input_file2"
			_polap_log2 "    concatenated fastq file: $input_file1 + $input_file2 -> $output_file"
		else
			_polap_log2 "    concatenated fastq file: $input_file1 -> $output_file"
		fi
	)

}

# return the type of pacbio fastq data type based on the average quality score
# input1: fastq
# output1: fastq.tsv
# output2: fastq.txt
_polap_lib_fastq-pacbio-type() {
	local input_file1="$1" # First input file
	local output_file="$2" # Output file (e.g., s1.fq)

	# Check if the first input file exists
	if [[ ! -f "$input_file1" ]]; then
		die "ERROR: no such file: $input_file1"
	fi

	if [[ ! -s "${input_file1}".tsv ]]; then
		seqkit head -n 10000 "${input_file1}" | seqkit stats -aT -o "${input_file1}.tsv"
	fi

	Rscript ${_POLAPLIB_DIR}/polap-r-fastq-pacbio-type.R \
		--tsv "${input_file1}.tsv" \
		-o "${output_file}"
	# -o "${input_file1}.txt" 2>/dev/null

}

# Version: v0.1.0
# Usage: seqkit_stats_outname /path/to/abc.fq[.gz] | /path/to/abc.fastq[.gz]
# Prints: /path/to/abc.fq.seqkit.stats.ta.txt
seqkit_stats_outname() {
	local infile="${1-}"
	if [[ -z "$infile" ]]; then
		echo "Usage: seqkit_stats_outname /path/to/abc.fq[.gz]|.fastq[.gz]" >&2
		return 2
	fi

	# Split path without external commands (safe under `set -euo pipefail`)
	local base="${infile##*/}" # filename
	local dir="${infile%/*}"   # directory (may be same as infile if no '/')
	[[ "$dir" == "$infile" ]] && dir="."

	# Normalize base name: drop .gz, then drop .fastq or .fq
	local stem="$base"
	stem="${stem%.gz}"
	stem="${stem%.fastq}"
	stem="${stem%.fq}"

	# Compose desired output name
	printf '%s\n' "${dir%/}/${stem}.fq.seqkit.stats.ta.txt"
}

# _polap_compress_or_decompress
# Version: v0.2.0
# Adaptive decompression / compression / linking with pv progress
# ---------------------------------------------------------------
# Behavior:
#   infile .gz  -> outfile plain  -> decompress (gzip -dc)
#   infile plain -> outfile .gz   -> compress   (gzip -c)
#   both .gz or both plain        -> make symlink (ln -sf)
#
# Example:
#   _polap_compress_or_decompress sample.fq.gz sample.fq
#   _polap_compress_or_decompress sample.fq sample.fq.gz
#   _polap_compress_or_decompress sample.fq.gz copy.fq.gz
#   _polap_compress_or_decompress sample.fq copy.fq
#
_polap_compress_or_decompress() {
	local infile="$1"
	local outfile="$2"

	[[ -r "$infile" ]] || {
		echo "[ERR] input not readable: $infile" >&2
		return 2
	}
	mkdir -p "$(dirname "$outfile")"

	# detect compression by extension
	local in_is_gz=0 out_is_gz=0
	[[ "$infile" =~ \.gz$ ]] && in_is_gz=1
	[[ "$outfile" =~ \.gz$ ]] && out_is_gz=1

	# pick pv (with file size) or cat
	local -a pv_cmd
	if command -v pv >/dev/null 2>&1; then
		local size
		size=$(stat -c %s "$infile" 2>/dev/null || stat -f %z "$infile" 2>/dev/null || echo "")
		if [[ -n "$size" ]]; then
			pv_cmd=(pv -pterb -s "$size" "$infile")
		else
			pv_cmd=(pv -pterb "$infile")
		fi
	else
		pv_cmd=(cat "$infile")
	fi

	_polap_log1 "[INFO] infile=$infile -> outfile=$outfile"

	# same compression type -> symlink
	if ((in_is_gz == out_is_gz)); then
		_polap_log2 "[INFO] same type -> creating symlink"
		ln -sf "$(realpath -s "$infile")" "$outfile"
		return 0
	fi

	# decompress
	if ((in_is_gz == 1 && out_is_gz == 0)); then
		_polap_log2 "[INFO] decompressing (.gz -> plain)"
		"${pv_cmd[@]}" | gzip -dc >"$outfile"
		return 0
	fi

	# compress
	if ((in_is_gz == 0 && out_is_gz == 1)); then
		_polap_log2 "[INFO] compressing (plain -> .gz)"
		"${pv_cmd[@]}" | gzip -c >"$outfile"
		return 0
	fi

	# fallback (weird combination)
	_polap_log2 "[WARN] unknown combination, copying"
	"${pv_cmd[@]}" >"$outfile"
}

# infile outfile size
_polap_lib_fastq-sample-to() {
	local infile="${1}"
	local outfile="${2}"
	local max_size="${3}"

	# check if seqkit stats -Ta output exists.
	# local infile_seqkit_stats=$(seqkit_stats_outname "$infile")
	# _polap_log2 "infile_seqkit_stats: $infile_seqkit_stats"
	# if [[ ! -s "$infile_seqkit_stats" ]]; then
	# 	seqkit stats -Ta "${infile}" -o "${infile_seqkit_stats}"
	# fi
	# local sum_size=$(cat "${infile_seqkit_stats}" | awk 'NR==2 {print $5}')
	local sum_size=$(_polap_lib_fastq-estimate-bases ${infile})

	_polap_log2 "sum_size: ${sum_size}"
	max_size=$(_polap_lib_unit-convert_to_int ${max_size})
	_polap_log2 "max_size: ${max_size}"

	local rate
	if [[ "${max_size}" == "0" ]]; then
		local rate=1
	else
		local rate=$(echo "scale=9; $max_size / $sum_size" | bc)
	fi
	_polap_log2 "rate: ${rate}"

	if [[ $(echo "$rate < 1" | bc) -eq 1 ]]; then

		_polap_lib_random-get

		local seed=${_polap_var_random_number}

		_polap_log2 "  input1: ${infile}"
		_polap_log2 "  output1: ${outfile}"
		_polap_log2 "  random seed: ${seed}"
		_polap_log2 "  long-read: ${max_size} (bp)"
		_polap_log2 "  sampling rate: ${rate}"

		_polap_log2 "sampling using seqkit sample with rate: ${rate}, max_size: ${max_size} and seed: ${seed} for ${infile} -> ${outfile}"

		if [[ -s "${infile}" ]]; then
			rm -f "${outfile}"
			seqkit sample \
				-p "$rate" \
				-s "${seed}" \
				"${infile}" \
				-o "${outfile}" \
				2>/dev/null
		else
			_polap_log0 "ERROR: no such file: ${infile}"
			exit 1
		fi
	else
		# differently for the output file name and the input file name
		_polap_compress_or_decompress "$infile" "$outfile"

		_polap_log2 "No sampling because the rate is greater than 1."
	fi
}

_polap_lib_fastq-sample() {
	local infile="${1}"
	local outfile="${2}"
	local rate="${3}"

	if [[ $(echo "$rate < 1" | bc) -eq 1 ]]; then
		_polap_lib_random-get
		local seed=${_polap_var_random_number}
		_polap_log2 "sampling using seqkit sample with rate: ${rate} and seed: ${seed} for ${infile} -> ${outfile}"
		if [[ -s "${infile}" ]]; then
			rm -f "${outfile}"
			seqkit sample \
				-p "$rate" \
				-s "${seed}" \
				"${infile}" \
				-o "${outfile}" \
				2>/dev/null
		else
			_polap_log0 "ERROR: no such file: ${infile}"
			exit 1
		fi
	else
		# ln -fs $(basename "${infile}") "${outfile}"
		# ln -fs "${infile}" "${outfile}"
		_polap_lib_filepath-smart_ln_s2 "${infile}" "${outfile}"
		_polap_log2 "No sampling because the rate is greater than 1: for ${infile} -> ${outfile}"
	fi
}

# Use the first 25 sequences to determine the sequencing data type
# illumina: average length < 500
# pacbio-hifi: length >= 500 and Q30 greater than 90.
# nano-raw: length >= 500 and Q30 greater than 30.
_polap_lib_fastq-check-type() {
	local file="$1"

	if [[ ! -f "$file" ]]; then
		# echo "Error: File not found: $file" >&3
		echo "unknown"
		return 1
	fi

	# Use head to avoid reading the full file
	local stats
	stats=$(head -n 20000 "$file" | tail -n 10000 | seqkit stats -Ta 2>/dev/null | awk 'NR==2')

	if [[ -z "$stats" ]]; then
		echo "Error: seqkit stats failed on $file" >&2
		echo "unknown"
		return 1
	fi

	local avg_len q30
	avg_len=$(echo "$stats" | awk '{print $7}')
	q30=$(echo "$stats" | awk '{print $15}')
	avgqual=$(echo "$stats" | awk '{print $16}')

	if (($(echo "$avg_len < 500" | bc -l))); then
		echo "illumina"
	elif (($(echo "$avg_len >= 500 && $q30 > 90 && $avgqual > 20" | bc -l))); then
		echo "pacbio-hifi"
	elif (($(echo "$avg_len >= 500 && $q30 < 50 && $avgqual < 20" | bc -l))); then
		echo "nano-raw"
	else
		echo "unknown"
	fi
}

_polap_lib_fastq-sample-to-coverage() {
	local infile="${1}"
	local outfile="${2}"
	local coverage="${3}"
	local genomesize="${4}"

	_polap_log1 "subsample the long-read data using a given target coverage: ${coverage}x"
	_polap_lib_fastq-total-length-of "${infile}" "${_arg_outdir}/l.fq.txt"
	local _l=$(<"${_arg_outdir}/l.fq.txt")
	local _v="${genomesize}"

	local _genome_coverage=$(echo "scale=5; ${_l} / ${_v}" | bc)
	local _rate=$(echo "scale=5; ${coverage} / ${_genome_coverage}" | bc)

	_polap_lib_random-get
	local _seed=${_polap_var_random_number}

	_polap_log1 "  input1: ${infile}"
	_polap_log1 "  output1: ${outfile}"
	_polap_log1 "  random seed: ${_seed}"
	_polap_log1 "  long-read: ${_l} (bp)"
	_polap_log1 "  genome size: ${_v} (bp)"
	_polap_log1 "  long-read genome coverage: ${_genome_coverage}x"
	_polap_log1 "  target coverage: ${coverage}x"
	_polap_log1 "  sampling rate: ${_rate}"

	local result=$(echo "$_rate < 1" | bc)

	if [ "$result" -eq 1 ]; then
		# echo "The rate value is less than 1"
		if [[ "${_arg_dry}" == "off" ]]; then
			rm -f "${outfile}"
			seqkit sample \
				-p "${_rate}" \
				-s "${_seed}" \
				"${infile}" \
				-o ${outfile} 2>${_polap_output_dest}
			# gzip "${_outfile}"
		fi

	else
		# echo "The value is not less than 1"
		_polap_log1 "  sampling rate is not less than 1: ${_rate}"
		_polap_log1 "  no subsampling of input: ${infile}"
		_polap_log0 "  no subsampling: ${outfile}"
		rm -f "${outfile}"
		# ln -s "$(realpath ${_infile})" "$(realpath -m ${_outfile})"
		_polap_lib_filepath-smart_ln_s2 "${infile}" "${outfile}"
		# _polap_lib_make_relative_symlink "${_infile}" "${_outfile}"
	fi

	return 0
}

_polap_lib_fastq-normalize-filename() {
	local infile="$1"
	local outdir="${_arg_outdir}/tmp"
	mkdir -p "$outdir"

	# Get base name without path
	local base="$(basename "$infile")"
	local stem=l

	case "$base" in
	*.fq | *.fastq)
		cp "$infile" "$outdir/${stem}.fq"
		;;
	*.fq.gz | *.fastq.gz)
		gunzip -c "$infile" >"$outdir/${stem}.fq"
		;;
	*.fq.tar.gz | *.fastq.tar.gz)
		# Extract into temp dir
		tmpdir=$(mktemp -d)
		tar -xzf "$infile" -C "$tmpdir"
		# Concatenate all fastq/fq files into one
		cat "$tmpdir"/*.{fq,fastq} 2>/dev/null >"$outdir/${stem}.fq"
		rm -rf "$tmpdir"
		;;
	*.fa | *.fasta)
		cp "$infile" "$outdir/${stem}.fa"
		;;
	*.fa.gz | *.fasta.gz)
		gunzip -c "$infile" >"$outdir/${stem}.fa"
		;;
	*.fa.tar.gz | *.fasta.tar.gz)
		tmpdir=$(mktemp -d)
		tar -xzf "$infile" -C "$tmpdir"
		cat "$tmpdir"/*.{fa,fasta} 2>/dev/null >"$outdir/${stem}.fa"
		rm -rf "$tmpdir"
		;;
	*)
		echo "Unsupported file type: $infile" >&2
		return 1
		;;
	esac
}

# Version: v0.2.0
# Estimate total bases in a FASTQ (compressed or not) from file size.
# Prints only the numeric total base count (integer).
# Usage:
#   _polap_lib_fastq-estimate-bases reads.fastq[.gz|.bz2|.xz]

_polap_lib_fastq-estimate-bases() {
	local fq="$1"
	local ratio=1.0 comp bytes bases

	# Get file size
	bytes=$(stat -L -c %s "$fq" 2>/dev/null) || {
		echo "0" >&2
		return 1
	}

	# Detect compression type
	case "$fq" in
	*.gz)
		comp="gz"
		ratio=3.5
		;; # default ONT
	*.bz2)
		comp="bz2"
		ratio=3.2
		;;
	*.xz)
		comp="xz"
		ratio=4.0
		;;
	*)
		comp="none"
		ratio=1.0
		;;
	esac

	# Estimate total bases = bytes × ratio × 0.45
	bases=$(echo "$bytes * $ratio * 0.45" | bc)

	# Print only the number (integer)
	printf "%.0f\n" "$bases"
}
