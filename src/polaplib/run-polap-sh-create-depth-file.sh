#!/usr/bin/env bash
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
# Bash script that prints 2 lines in a TSV file.
# It is used in the seed contig selection procedure.
# This would create the output file with the depth range given by the first and
# the second argumnt. The third argument is the output text filename.
#
# Example:
# <cmd> 1 2 out.tsv
# Output:
# depth_lower_bound	depth_upper_bound
# 1	2
#
# See Also:
# polap-function-disassemble-seeds.sh
# run-polap-function-seeds.sh
#
# TODO: rename: polap-bash-create-depth-file.sh
#
# Check: 2025-06-17
################################################################################

# Only 3 arguments
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 lower_bound upper_bound output_file"
  exit 1
fi

lower_bound=$1
upper_bound=$2
output_file=$3

# If -e is in effect, the following sequences are recognized:
# \n     new line
# \r     carriage return
# \t     horizontal tab
# In MacOS, we do not need -e option, but in Linux, we need -e option to use
# the escape character.
# In MacOS, the -e option is ignored. So, this code works for MacOS and Linux.
echo -e "depth_lower_bound\tdepth_upper_bound" >"$output_file"
echo -e "$lower_bound\t$upper_bound" >>"$output_file"

