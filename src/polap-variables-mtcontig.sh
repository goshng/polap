# base folders
local _polap_var_ga="${ODIR}/${INUM}"
local _polap_var_mtcontigs="${_polap_var_ga}/${JNUM}"/mtcontigs

# input
local _polap_var_assembly_graph_final_gfa="${_polap_var_ga}/30-contigger/graph_final.gfa"
local _polap_var_annotation_table="${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt"

# output
local MTCONTIGNAME="${_polap_var_ga}"/mt.contig.name-"${JNUM}"
local _polap_var_mtcontig_table="${_polap_var_mtcontigs}/1-mtcontig.table.tsv"

# run-polap-select-contigs-by-1-annotation.R
local _polap_var_mtcontig_base="${_polap_var_mtcontigs}/1-mtcontig"
local _polap_var_mtcontig_annotated="${_polap_var_mtcontigs}/1-mtcontig.annotated.txt"
local _polap_var_mtcontig_annotated_stats="${_polap_var_mtcontigs}/1-mtcontig.annotated.stats.txt"
local _polap_var_mtcontig_mixfit="${_polap_var_mtcontigs}/1-mtcontig.mixfit.txt"

local _polap_var_mtcontig_depth_stats="${_polap_var_mtcontigs}/1-mtcontig.depth.stats.txt"

# local _polap_var_mtcontigs_mt_stats="${_polap_var_mtcontigs}/1-mtcontig.mt.stats.txt"
# local _polap_var_mtcontigs_pt_stats="${_polap_var_mtcontigs}/1-mtcontig.pt.stats.txt"

# run-polap-select-contigs-by-3-subset-by-depth.R
local _polap_var_gfa_all="${_polap_var_mtcontigs}/2-gfa.all.gfa"
local _polap_var_gfa_seq_part="${_polap_var_mtcontigs}/2-gfa.seq.part.tsv"
local _polap_var_gfa_seq_filtered="${_polap_var_mtcontigs}/2-gfa.seq.filtered.txt"
local _polap_var_gfa_seq_filtered_range="${_polap_var_mtcontigs}/2-gfa.seq.filtered.range.txt"

#
local _polap_var_gfa_seq_filtered_edge="${_polap_var_mtcontigs}/2-gfa.seq.filtered.edge.txt"
local _polap_var_gfa_filtered="${_polap_var_mtcontigs}/2-gfa.filtered.gfa"
local _polap_var_gfa_links="${_polap_var_mtcontigs}/3-gfa.links.tsv"
local _polap_var_gfa_links_number="${_polap_var_mtcontigs}/3-gfa.links.number.txt"
local _polap_var_gfa_links_order="${_polap_var_mtcontigs}/3-gfa.links.order.txt"
local _polap_var_gfa_links_contig="${_polap_var_mtcontigs}/3-gfa.links.contig.txt"
local _polap_var_gfa_links_contig_na="${_polap_var_mtcontigs}/3-gfa.links.contig.na.txt"
local _polap_var_gfa_links_seed="${_polap_var_mtcontigs}/4-gfa.links.seed.txt"
local _polap_var_gfa_links_mtcontig="${_polap_var_mtcontigs}/5-gfa.links.mtcontig.txt"
