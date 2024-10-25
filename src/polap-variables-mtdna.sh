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

# File paths
local _polap_var_assembly_graph_gfa="${_polap_var_ga}/assembly_graph.gfa"
local _polap_var_annotation_table="${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt"
local _polap_var_mt_fasta="${_polap_var_ga}/mt.0.fasta"
local _polap_var_mt_edges="${_polap_var_ga}/mt.0.edges"

# _polap_var_mtdna_
local _polap_var_1_gfa_all="${_polap_var_mtdna}/1-gfa.all.gfa"
local _polap_var_gfa_links="${_polap_var_mtdna}/1-gfa.links.tsv"
local _polap_var_gfa_links_edges="${_polap_var_mtdna}/1-gfa.links.edges.txt"
local _polap_var_gfa_links_circular_path="${_polap_var_mtdna}/2-gfa.links.circular.path.txt"
local _polap_var_circular_path="${_polap_var_mtdna}/3-circular.path.txt"
local _polap_var_edge_fasta="${_polap_var_mtdna}/5-edge.fasta"
local _polap_var_fasta="${_polap_var_mtdna}/4-gfa.fasta"
