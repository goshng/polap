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
# The following two scripts used to be used for archiving output files.
# Now, we use rsync and archive template text file for archiving.
# We could delete them later but leave them now.
# polaplib/polap-package-common.sh
# polaplib/polap-package-mtcontigs.sh
################################################################################

# common: mtcontigs
local _ppack_var_ga_mtcontigs="${_ppack_var_ga}/51-mtcontigs"
local _ppack_var_mtcontigs="${_ppack_var_ga}/51-mtcontigs/${_arg_knum}"
local _ppack_var_mtcontigs_depth_range_preselection="${_ppack_var_mtcontigs}/depth.range.preselection.txt"
local _ppack_var_mtcontigs_depth_range_graphfilter="${_ppack_var_mtcontigs}/depth.range.graphfilter.txt"
local _ppack_var_mtcontigs_1_custom_depth_range="${_ppack_var_mtcontigs}/custom.depth.range-1.txt"
local _ppack_var_mtcontigs_2_custom_depth_range="${_ppack_var_mtcontigs}/custom.depth.range-2.txt"
local _ppack_var_mtcontigs_preselection="${_ppack_var_mtcontigs}/1-preselection.by.gene.density.txt"
local _ppack_var_mtcontigs_gfa_all="${_ppack_var_mtcontigs}/3-gfa.all.gfa"
local _ppack_var_mtcontigs_gfa_seq_filtered="${_ppack_var_mtcontigs}/3-gfa.seq.depthfiltered.txt"
local _ppack_var_mtcontigs_gfa_seq_part="${_ppack_var_mtcontigs}/3-gfa.seq.all.tsv"
local _ppack_var_mtcontigs_gfa_seq_filtered_edge="${_ppack_var_mtcontigs}/4-gfa.seq.depthfiltered.edge.txt"
local _ppack_var_mtcontigs_gfa_depthfiltered_gfa="${_ppack_var_mtcontigs}/4-gfa.depthfiltered.gfa"

local _ppack_var_mtcontigs_links="${_ppack_var_mtcontigs}/4-gfa.links"
local _ppack_var_mtcontigs_links_tsv="${_ppack_var_mtcontigs}/4-gfa.links.tsv"
local _ppack_var_mtcontigs_links_contig_na="${_ppack_var_mtcontigs}/4-gfa.links.contig.na.txt"
local _ppack_var_mtcontigs_links_contig="${_ppack_var_mtcontigs}/4-gfa.links.contig.txt"
local _ppack_var_mtcontigs_links_number="${_ppack_var_mtcontigs}/4-gfa.links.number.txt"
local _ppack_var_mtcontigs_links_order="${_ppack_var_mtcontigs}/4-gfa.links.order.txt"
local _ppack_var_mtcontigs_links_seed="${_ppack_var_mtcontigs}/5-gfa.links.seed.txt"
local _ppack_var_mtcontigs_links_mtcontig="${_ppack_var_mtcontigs}/6-gfa.links.mtcontig.txt"
local _ppack_var_mtcontigs_7mtcontigname="${_ppack_var_mtcontigs}/7-mt.contig.name.txt"
local _ppack_var_mtcontigs_8mtcontigname="${_ppack_var_mtcontigs}/8-mt.contig.name.txt"
local _ppack_var_mtcontig_table="${_ppack_var_mtcontigs}/8-mtcontig.table.tsv"
local _ppack_var_mtcontigs_annotation_table_seed="${_ppack_var_mtcontigs}/8-mtcontig-annotation-table-seed.txt"
# delete these later
local _ppack_var_mtcontigs_2_depth_range_by_cdf_copy_number="${_ppack_var_mtcontigs}/2-depth.range.by.cdf.copy.number.txt"
local _ppack_var_mtcontig_table="${_ppack_var_mtcontigs}/8-mtcontig.table.tsv"
# delete these later
local _ppack_var_mtcontigs_2_depth_range_by_cdf_copy_number="${_ppack_var_mtcontigs}/2-depth.range.by.cdf.copy.number.txt"
