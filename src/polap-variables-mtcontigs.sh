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

# common: mtcontigs
local _polap_var_ga_mtcontigs="${_polap_var_ga}/51-mtcontigs"
local _polap_var_mtcontigs="${_polap_var_ga}/51-mtcontigs/${KNUM}"
local _polap_var_mtcontigs_depth_range_preselection="${_polap_var_mtcontigs}/depth.range.preselection.txt"
local _polap_var_mtcontigs_depth_range_graphfilter="${_polap_var_mtcontigs}/depth.range.graphfilter.txt"
local _polap_var_mtcontigs_1_custom_depth_range="${_polap_var_mtcontigs}/custom.depth.range-1.txt"
local _polap_var_mtcontigs_2_custom_depth_range="${_polap_var_mtcontigs}/custom.depth.range-2.txt"
local _polap_var_mtcontigs_preselection="${_polap_var_mtcontigs}/1-preselection.by.gene.density.txt"
local _polap_var_mtcontigs_gfa_all="${_polap_var_mtcontigs}/3-gfa.all.gfa"
local _polap_var_mtcontigs_gfa_seq_filtered="${_polap_var_mtcontigs}/3-gfa.seq.depthfiltered.txt"
local _polap_var_mtcontigs_gfa_seq_part="${_polap_var_mtcontigs}/3-gfa.seq.all.tsv"
local _polap_var_mtcontigs_gfa_seq_filtered_edge="${_polap_var_mtcontigs}/4-gfa.seq.depthfiltered.edge.txt"
local _polap_var_mtcontigs_gfa_depthfiltered_gfa="${_polap_var_mtcontigs}/4-gfa.depthfiltered.gfa"

local _polap_var_mtcontigs_links="${_polap_var_mtcontigs}/4-gfa.links"
local _polap_var_mtcontigs_links_tsv="${_polap_var_mtcontigs}/4-gfa.links.tsv"
local _polap_var_mtcontigs_links_contig_na="${_polap_var_mtcontigs}/4-gfa.links.contig.na.txt"
local _polap_var_mtcontigs_links_contig="${_polap_var_mtcontigs}/4-gfa.links.contig.txt"
local _polap_var_mtcontigs_links_number="${_polap_var_mtcontigs}/4-gfa.links.number.txt"
local _polap_var_mtcontigs_links_order="${_polap_var_mtcontigs}/4-gfa.links.order.txt"
local _polap_var_mtcontigs_links_seed="${_polap_var_mtcontigs}/5-gfa.links.seed.txt"
local _polap_var_mtcontigs_links_mtcontig="${_polap_var_mtcontigs}/6-gfa.links.mtcontig.txt"
local _polap_var_mtcontigs_7mtcontigname="${_polap_var_mtcontigs}/7-mt.contig.name.txt"
local _polap_var_mtcontigs_8mtcontigname="${_polap_var_mtcontigs}/8-mt.contig.name.txt"
local _polap_var_mtcontig_table="${_polap_var_mtcontigs}/8-mtcontig.table.tsv"
# delete these later
local _polap_var_mtcontigs_2_depth_range_by_cdf_copy_number="${_polap_var_mtcontigs}/2-depth.range.by.cdf.copy.number.txt"
local _polap_var_mtcontig_table="${_polap_var_mtcontigs}/8-mtcontig.table.tsv"
# delete these later
local _polap_var_mtcontigs_2_depth_range_by_cdf_copy_number="${_polap_var_mtcontigs}/2-depth.range.by.cdf.copy.number.txt"
