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
# This script is used for converting minimap2's PAF to TAB file.
#
# See Also:
# This script was used by map-read subcommand in run-polap-function-oga.
#
# TODO: rename: polap-bash-minimap2-paf2tab.sh
#
# Check: 2025-06-17
################################################################################

# only 3 arguments
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 min_read_length paf_file tab_file"
  exit 1
fi

_min_read_length=$1
_paf=$2
_tab=$3

# This filters PAF by the lengths of reads.
# We use only the first 11 columns by excluding the last one.
cut -f1-11 "${_paf}" |
  awk -v minlength="${_min_read_length}" \
    '{if ($2>=minlength) {print}}' \
    >"${_tab}"
