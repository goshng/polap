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
# Bash functions for file paths
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

# Portable function to resolve full script path including symlinks
_polap_lib_filepath-get_script_path() {
	local target=$1
	local prev_dir
	local dir

	# Resolve $target until it's not a symlink
	while [ -h "$target" ]; do
		dir="$(cd "$(dirname "$target")" && pwd)" || exit 1
		target="$(readlink "$target")"

		# If readlink gives a relative path, prepend the symlink's directory
		case "$target" in
		/*) ;; # Absolute path, leave it
		*) target="$dir/$target" ;;
		esac
	done

	# Final absolute path
	dir="$(cd "$(dirname "$target")" && pwd)" || exit 1
	echo "$dir/$(basename "$target")"
}

# Determine script path only if executed (not sourced)
# if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
# 	script_path="$(_polap_lib_filepath-get_script_path "$0")"
# 	script_dir="$(dirname "$script_path")"
#
# 	echo "Resolved script path: $script_path"
# 	echo "Script directory:    $script_dir"
# fi

_polap_lib_make_relative_symlink() {
	local target="$1"
	local linkpath="$2"

	if [[ -z "$target" || -z "$linkpath" ]]; then
		echo "Usage: make_relative_symlink <target-file> <symlink-path>" >&2
		return 1
	fi

	# Resolve absolute target path
	local abs_target
	abs_target=$(realpath "$target") || {
		echo "Error: target '$target' does not exist." >&2
		return 2
	}

	# Create parent directory for link if needed
	local link_dir
	link_dir=$(dirname "$linkpath")
	mkdir -p "$link_dir" || {
		echo "Error: Failed to create parent directory '$link_dir'" >&2
		return 3
	}

	# Compute relative path from link location to target
	local rel_target
	rel_target=$(realpath --relative-to="$link_dir" "$abs_target") || {
		echo "Error: Failed to compute relative path." >&2
		return 4
	}

	# Remove existing file/symlink if exists
	if [[ -e "$linkpath" || -L "$linkpath" ]]; then
		rm -f "$linkpath" || {
			echo "Error: Failed to remove existing '$linkpath'" >&2
			return 5
		}
	fi

	# Create the symlink
	ln -s "$rel_target" "$linkpath" || {
		echo "Error: Failed to create symlink '$linkpath'" >&2
		return 6
	}

	echo "Symlink created: $linkpath -> $rel_target"
}

_polap_lib_filepath-smart_ln_s() {
	local src="$1"
	local tgt="$2"
	local src_abs tgt_abs src_dir tgt_dir

	src_abs=$(realpath -m "$src")
	tgt_abs=$(realpath -m "$tgt")

	src_dir=$(dirname "$src_abs")
	tgt_dir=$(dirname "$tgt_abs")

	# When source and target are in the same directory
	if [[ "$src_dir" == "$tgt_dir" ]]; then
		ln -fs "$(basename "$src")" "$tgt"
	else
		# Try to compute relative path from target to source
		rel_src=$(realpath --relative-to="$tgt_dir" "$src_abs")
		ln -fs "$rel_src" "$tgt"
	fi
}

_polap_lib_filepath-smart_ln_s2() {
	local src="$1"
	local tgt="$2"
	local src_abs src_dir tgt_dir rel_src

	src_abs=$(realpath -m "$src") || return 1
	src_dir=$(dirname "$src_abs")

	tgt_dir=$(realpath -m "$(dirname "$tgt")") || return 1

	# When source and target are in the same directory
	if [[ "$src_dir" == "$tgt_dir" ]]; then
		ln -fs "$(basename "$src")" "$tgt"
	else
		# Compute relative path from target dir to source
		rel_src=$(realpath --relative-to="$tgt_dir" "$src_abs")
		ln -fs "$rel_src" "$tgt"
	fi
}
