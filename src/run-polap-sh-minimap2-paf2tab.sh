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

# Check if exactly 3 arguments are provided
if [ "$#" -ne 3 ]; then
	echo "Usage: $0 min_read_length paf_file tab_file"
	exit 1
fi

# Assign positional arguments to variables
_min_read_length=$1
_paf=$2
_tab=$3

# Create the output file with the depth range
cut -f1-11 "${_paf}" |
	awk -v minlength="${_min_read_length}" \
	'{if ($2>=minlength) {print}}' \
		>"${_tab}"
