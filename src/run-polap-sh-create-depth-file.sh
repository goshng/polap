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
	echo "Usage: $0 lower_bound upper_bound output_file"
	exit 1
fi

# Assign positional arguments to variables
lower_bound=$1
upper_bound=$2
output_file=$3

# Create the output file with the depth range
echo -e "depth_lower_bound\tdepth_upper_bound" >"$output_file"
echo -e "$lower_bound\t$upper_bound" >>"$output_file"
