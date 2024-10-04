# Set variables
local _polap_var_oga="${ODIR}/${INUM}"
local _polap_var_mtdna="${_polap_var_oga}/mtdna"

# File paths
local _polap_var_assembly_graph_final_gfa="${_polap_var_oga}/assembly_graph.gfa"
local _polap_var_annotation_table="${_polap_var_oga}/assembly_info_organelle_annotation_count-all.txt"
local _polap_var_mt_fasta="${_polap_var_oga}/mt.0.fasta"
local _polap_var_mt_edges="${_polap_var_oga}/mt.0.edges"
local _polap_var_gfa_all="${_polap_var_mtdna}/1-gfa.all.gfa"

local _polap_var_gfa_links="${_polap_var_mtdna}/1-gfa.links.tsv"
local _polap_var_gfa_links_edges="${_polap_var_mtdna}/1-gfa.links.edges.txt"
local _polap_var_gfa_links_circular_path="${_polap_var_mtdna}/2-gfa.links.circular.path.txt"
local _polap_var_circular_path="${_polap_var_mtdna}/3-circular.path.txt"
local _polap_var_edge_fasta="${_polap_var_mtdna}/5-edge.fasta"
local _polap_var_fasta="${_polap_var_mtdna}/4-gfa.fasta"
