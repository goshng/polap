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

# base folders
local _polap_var_ga="${ODIR}/${INUM}"
local _polap_var_manual_depth_range="${ODIR}"/0/1-manual.depth.range.txt
local _polap_var_manual_copy_range="${ODIR}"/0/1-manual.copy.range.txt
local _polap_var_mtcontigs="${_polap_var_ga}/${JNUM}"/mtcontigs

local _polap_var_1_custom_depth_range="${ODIR}"/0/1-custom.depth.range.txt
local _polap_var_2_custom_depth_range="${ODIR}"/0/2-custom.depth.range.txt

local _polap_var_2_depth_range_by_cdf_copy_number="${ODIR}/0/2-depth.range.by.cdf.copy.number.txt"
local _polap_var_3_depth_range_by_mixture="${ODIR}/0/3-depth.range.by.mixture.txt"
local _polap_var_3_mixfit="${_polap_var_ga}/3-mixfit.txt"

# input
local _polap_var_assembly_graph_final_gfa="${_polap_var_ga}/30-contigger/graph_final.gfa"
local _polap_var_annotation_table="${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt"

# output
local MTCONTIGNAME="${_polap_var_ga}"/mt.contig.name-"${JNUM}"
local _polap_var_mtcontig_table="${_polap_var_mtcontigs}/1-mtcontig.table.tsv"

# run-polap-r-select-contigs-by-1-preselection.R
local _polap_var_preselection_by_gene_density="${_polap_var_mtcontigs}/1-preselection.by.gene.density.txt"
local _polap_var_preselection_by_depth_mixture="${_polap_var_mtcontigs}/1-preselection.by.depth.mixture.txt"
local _polap_var_depth_range_by_gene_density="${_polap_var_mtcontigs}/1-depth.range.by.gene.density.txt"
local _polap_var_depth_range_by_depth_mixture="${_polap_var_mtcontigs}/1-depth.range.by.depth.mixture.txt"
local _polap_var_mixfit="${_polap_var_mtcontigs}/1-mixfit.txt"

# run-polap-r-select-contigs-by-2-determine-depth-range.R
local _polap_var_depth_range_by_cdf_copy_number="${_polap_var_mtcontigs}/2-depth.range.by.cdf.copy.number.txt"

# run-polap-select-contigs-by-1-annotation.R
local _polap_var_mtcontig_depth_range="${_polap_var_mtcontigs}/1-depth.range.txt"

local _polap_var_mtcontig_base="${_polap_var_mtcontigs}/1-mtcontig"
local _polap_var_mtcontig_annotated="${_polap_var_mtcontigs}/1-mtcontig.annotated.txt"
local _polap_var_mtcontig_annotated_stats="${_polap_var_mtcontigs}/1-mtcontig.annotated.stats.txt"

# local _polap_var_mtcontigs_mt_stats="${_polap_var_mtcontigs}/1-mtcontig.mt.stats.txt"
# local _polap_var_mtcontigs_pt_stats="${_polap_var_mtcontigs}/1-mtcontig.pt.stats.txt"

# Step 3: filtering GFA using the depth range from step 2
# run-polap-select-contigs-by-3-subset-by-depth.R
local _polap_var_gfa_all="${_polap_var_mtcontigs}/3-gfa.all.gfa"
local _polap_var_gfa_seq_part="${_polap_var_mtcontigs}/3-gfa.seq.part.tsv"
local _polap_var_gfa_seq_filtered="${_polap_var_mtcontigs}/3-gfa.seq.filtered.txt"
# local _polap_var_gfa_seq_filtered_range="${_polap_var_mtcontigs}/3-gfa.seq.filtered.range.txt"

# Step 4: preparing the graph for finding connected components
local _polap_var_gfa_seq_filtered_edge="${_polap_var_mtcontigs}/2-gfa.seq.filtered.edge.txt"
local _polap_var_gfa_filtered="${_polap_var_mtcontigs}/2-gfa.filtered.gfa"

local _polap_var_links="${_polap_var_mtcontigs}/4-gfa.links"
local _polap_var_links_tsv="${_polap_var_links}.tsv"
local _polap_var_links_number="${_polap_var_links}.number.txt"
local _polap_var_links_order="${_polap_var_links}.order.txt"
local _polap_var_links_contig="${_polap_var_links}.contig.txt"
local _polap_var_links_contig_na="${_polap_var_links}.contig.na.txt"

# Step 5: Find connected components using Python script
local _polap_var_links_seed="${_polap_var_mtcontigs}/5-gfa.links.seed.txt"

# Step 6: converting the depth-filtered contigs in edge with numbers
local _polap_var_links_mtcontig="${_polap_var_mtcontigs}/6-gfa.links.mtcontig.txt"
