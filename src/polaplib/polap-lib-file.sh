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
# Bash functions for handling regular files.
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

# arg1: source folder
# arg2: destination folder
# arg3: maximum size of files for rsync or copying
_polap_lib_file-rsync() {
	local input1="$1"
	local output="$2"
	local _max_size="${3:-3M}"

	input1="${input1%/}"
	output="${output%/}"
	if [[ -d "${output}" ]]; then
		_polap_log0 "ERROR: folder already exists: ${output}"
	else
		rsync -aq --max-size=${_max_size} "${input1}/" "${output}/"
	fi
}

# Example usage:
# archive-folder "Afolder" "Bfolder" "template.txt"
_polap_lib_file-archive-folder() {
	local src_dir="$1"
	local dest_dir="$2"
	local template_file="$3"

	if [[ ! -d "$src_dir" ]]; then
		_polap_log0 "Error: Source directory '$src_dir' does not exist."
		return 1
	fi

	if [[ ! -f "$template_file" ]]; then
		_polap_log0 "Error: Template file '$template_file' does not exist."
		return 1
	fi

	mkdir -p "$dest_dir"

	while IFS= read -r file; do
		if [[ "$file" == *"#"* ]]; then
			continue
		fi
		src_path="$src_dir/$file"
		dest_path="$dest_dir/$file"
		dest_folder="$(dirname "$dest_path")"

		if [[ -f "$src_path" ]]; then
			mkdir -p "$dest_folder"
			cp -p "$src_path" "$dest_path"
		else
			echo "Warning: File '$src_path' does not exist, skipping."
		fi
	done <"$template_file"

	# Create a tar.gz archive of the destination directory
	tar -czf "${dest_dir}.tar.gz" -C "$(dirname "$dest_dir")" "$(basename "$dest_dir")"

	_polap_log1 "Archive created: ${dest_dir}.tar.gz"

	# Delete the destination directory after archiving
	# rm -rf "$dest_dir"
	# echo "Deleted temporary directory: $dest_dir"
}

# if _polap_lib_file-is_at_least_1MB "${_brg_outdir_i}/ont.asm.ec.fq"; then
#   echo "found: ${_brg_outdir_i}/ont.asm.ec.fq"
# else
_polap_lib_file-is_at_least_1MB() {
	local file="$1"
	[[ -f "$file" && $(stat -c%s "$file") -ge 1048576 ]]
}

# from disassemble.sh
#
# unzip a gzipped file leaving the input as is.
# do not unzip if the input is not a gzipped file.
#
# unzipped_file=$(_polap_gunzip_file "${_short_read1}")
# _rstatus="$?"
# if [[ "$_rstatus" -eq 0 ]]; then
#   _polap_log2 "  unzipped file: $unzipped_file"
#   _short_read1="$unzipped_file"
# fi
_polap_lib_file-gunzip() {
	local input_file="$1"

	# Check if the input file exists
	if [[ ! -f "$input_file" ]]; then
		die "Error: File '$input_file' not found."
	fi

	# Check if the file is gzipped
	if file "$input_file" | grep -q "gzip compressed data"; then
		# Extract the file name without the .gz extension
		local output_file="${input_file%.gz}"

		# Unzip the file and keep the original
		if gunzip -c "$input_file" >"$output_file"; then
			echo "$output_file"
			return 0
		else
			die "ERROR: failed to unzip '$input_file'."
		fi
	else
		_polap_log3 "    file '$input_file' is not gzipped."
		return 1
	fi
}
