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

local _polap_var_ga="${ODIR}/${INUM}"
local MTCONTIGNAME="${_polap_var_ga}"/mt.contig.name-"${JNUM}"
local _polap_var_ga_annotation="${_polap_var_ga}/assembly_info_organelle_annotation_count.txt"
local _polap_var_ga_annotation_all="${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt"
local _polap_var_ga_annotation_all_backup="${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt.backup"
local _polap_var_ga_annotation_table="${_polap_var_ga}/contig-annotation-table.txt"
local _polap_var_ga_annotation_depth_table="${_polap_var_ga}/contig-annotation-depth-table.txt"
local _polap_var_contigger="${_polap_var_ga}/30-contigger"
local _polap_var_contigger_gfa="${_polap_var_contigger}/graph_final.gfa"
local _polap_var_contigger_edges_stats="${_polap_var_contigger}/edges_stats.txt"
local _polap_var_contigger_edges_fasta="${_polap_var_contigger}/graph_final.fasta"
local _polap_var_contigger_contigs_stats="${_polap_var_contigger}/contigs_stats.txt"
local _polap_var_contigger_contigs_fasta="${_polap_var_contigger}/contigs.fasta"
local _polap_var_assembly_graph_final_gfa="${_polap_var_ga}/30-contigger/graph_final.gfa"
# temporary files to create edges_stats.txt
local _polap_var_ga_mtcontigs="${_polap_var_ga}/mtcontigs"
local _polap_var_ga_gfa_all="${_polap_var_ga_mtcontigs}/3-gfa.all.gfa"
local _polap_var_ga_gfa_seq_part="${_polap_var_ga_mtcontigs}/3-gfa.seq.part.tsv"
# assembly_contigs_stats="$FDIR"/30-contigger/contigs_stats.txt
# assembly_contigs_fasta="$FDIR"/30-contigger/contigs.fasta
local _polap_var_ann="${_polap_var_ga}"/50-annotation
local _polap_var_ann_CONTIGFILE="${_polap_var_ann}"/contig.fasta
local _polap_var_ann_CONTIGDB="${_polap_var_ann}"/contig
local _polap_var_ann_CONTIGNAME="${_polap_var_ann}"/contig.name
# MT
local _polap_var_ann_MTAABLAST="${_polap_var_ann}"/mtaa.blast
local _polap_var_ann_MTAABLASTBED="${_polap_var_ann}"/mtaa.blast.bed
local _polap_var_ann_MTAABED="${_polap_var_ann}"/mtaa.bed
local _polap_var_ann_MTGENECOUNT="${_polap_var_ann}"/mt.gene.count
# PT
local _polap_var_ann_PTAABLAST="${_polap_var_ann}"/ptaa.blast
local _polap_var_ann_PTAABLASTBED="${_polap_var_ann}"/ptaa.blast.bed
local _polap_var_ann_PTAABED="${_polap_var_ann}"/ptaa.bed
local _polap_var_ann_PTGENECOUNT="${_polap_var_ann}"/pt.gene.count
