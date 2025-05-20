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
