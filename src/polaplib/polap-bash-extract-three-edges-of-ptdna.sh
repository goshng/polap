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
# This script extracts the three edges from the plastid genome assembly graph.
# It is used in polap-cmd-disassemble.sh to extract them from ptGAUL
# assembly graph.
#
# Example:
# <cmd> assembly_graph.gfa mt.contig.name-1
# input in the gfa:
# P	contig_5	edge_4-,edge_3-,edge_5+,edge_3+,edge_4+,edge_3-	*
# output:
# edge_3
# edge_4
# edge_5
#
# Tips:
# echo edge_4-,edge_3-,edge_5+,edge_3+,edge_4+,edge_3- | tr ',' '\n' | sed 's/[+-]//' | sort -u
#
# See Also:
# polap-cmd-disassemble.sh
#
# Look for lines:
# extract three unique edges from the gfa for a ptDNA candidate
#
# L	edge_3	+	edge_4	+	0M	RC:i:41
# L	edge_3	+	edge_4	-	0M	RC:i:40
# L	edge_3	-	edge_5	-	0M	RC:i:42
# L	edge_3	-	edge_5	+	0M	RC:i:41
# P	contig_5	edge_4-,edge_3-,edge_5+,edge_3+,edge_4+,edge_3-	*
# P	contig_1	edge_1+	*
# P	contig_2	edge_2+	*
#
# We use only the contig with 3 unique edges. So, contig_1 and contig_2 are
# not used because they have only one edge.
# If a contig has exactly 3 unique edges, we print the 3 edges to the output
# file.
#
# TODO: rename: polap-bash-extract-three-edges-of-ptdna.sh
#
# Check: 2025-06-17
################################################################################

extract_edges() {
  local input_file="$1"
  local output_file="$2"

  # The output file starts with an empty content.
  >"$output_file"

  # Look for the contig lines with a starting P.
  ## P	contig_5	edge_4-,edge_3-,edge_5+,edge_3+,edge_4+,edge_3-	*
  # Read the contig and its following edges.
  ## contig: contig_5
  ## edges: edge_4-,edge_3-,edge_5+,edge_3+,edge_4+,edge_3-
  awk '$1 == "P" {print $2, $3}' "$input_file" | while read -r contig edges; do
    # edge numbers without signs and get the
    unique_edges=$(echo "$edges" | tr ',' '\n' | sed 's/[+-]//' | sort -u)

    edge_count=$(echo "$unique_edges" | wc -l)

    if [ "$edge_count" -eq 3 ]; then
      echo "$unique_edges" >>"$output_file"
    fi
  done
}

# Only 2 arguments
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 assembly_graph.gfa mt.contig.name-1"
  exit 1
fi

extract_edges "$1" "$2"
