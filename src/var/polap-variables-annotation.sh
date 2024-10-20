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

source "$script_dir/polap-variables-common.sh"

# local _polap_var_wga="${ODIR}/0"
local _polap_var_wga_contigger="${_polap_var_wga}/30-contigger"
local _polap_var_wga_contigger_gfa="${_polap_var_wga_contigger}/graph_final.gfa"
local _polap_var_wga_annotation="${_polap_var_wga}/assembly_info_organelle_annotation_count-all.txt"

ANUM=$INUM

local _polap_var_ga="$ODIR"/$INUM

assembly_contigs_stats="$FDIR"/30-contigger/contigs_stats.txt
assembly_contigs_fasta="$FDIR"/30-contigger/contigs.fasta

local _polap_var_ann="${_polap_var_ga}"/50-annotation
local _polap_var_ann_CONTIGFILE="${_polap_var_ann}"/contig.fasta
local _polap_var_ann_CONTIGDB="${_polap_var_ann}"/contig
local _polap_var_ann_MTAABLAST="${_polap_var_ann}"/mtaa.blast
local _polap_var_ann_MTAABED="${_polap_var_ann}"/mtaa.bed
local _polap_var_ann_MTGENECOUNT="${_polap_var_ann}"/mt.gene.count
local _polap_var_ann_PTAABLAST="${_polap_var_ann}"/ptaa.blast
local _polap_var_ann_PTAABED="${_polap_var_ann}"/ptaa.bed
local _polap_var_ann_PTGENECOUNT="${_polap_var_ann}"/pt.gene.count
local _polap_var_ann_CONTIGNAME="${_polap_var_ann}"/contig.name
