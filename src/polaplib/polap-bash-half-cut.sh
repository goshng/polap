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
# This script is used for circularizing a single contig.
# Because polap's contig mapping needs at least two seed edges to select reads,
# we cut a single contig into 2 pieces before mapping and selecting reads.
#
# See Also:
# This script was used in run-polap-function-oga.
#
# TODO: rename: polap-bash-half-cut.sh
#
# Check: 2025-06-17
################################################################################

# only 3 arguments
if [ "$#" -ne 3 ]; then
	echo "Usage: $0 mtdir mtdir/contig.fa mt.contig.name-1"
  echo "  output1: mt.contig.name-1-backup"
  echo "  output2: (new) mt.contig.name-1"
  echo "  output3: mtdir/contig.fa-backup"
  echo "  output4: (new) mtdir/contig.fa"
	exit 1
fi

_mtdir=$1
_mtdir_contig_fa=$2
_mt_contig_name=$3

# The following code used to be in run-polap-function-oga.sh. To simplify the code in there,
# we have this separate script that we execute in run-polap-function-oga.sh.
# We use it is subcommand map-reads.
cp "${_mtdir_contig_fa}" "${_mtdir_contig_fa}-backup"
seqkit fx2tab --length --name "${_mtdir_contig_fa}" -o "${_mtdir}"/contig.fa.len >/dev/null 2>&1
A=$(cut -f2 "${_mtdir}"/contig.fa.len)
B=$(echo "scale=0; $A/2" | bc)
C=$((B + 1))
seqkit subseq -r 1:"$B" "${_mtdir_contig_fa}" -o "${_mtdir}"/c1.fa >/dev/null 2>&1
seqkit subseq -r "$C":"$A" "${_mtdir_contig_fa}" -o "${_mtdir}"/c2.fa >/dev/null 2>&1
cat "${_mtdir}"/c?.fa | seqkit replace -p '.+' -r 'edge_{nr}' -o "${_mtdir_contig_fa} >/dev/null 2>&1
cp "${_mt_contig_name}" "${_mt_contig_name}-backup"
echo -e "edge_1\nedge_2" >"${_mt_contig_name}"
