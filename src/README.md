# 2024-10-20

Learn something about Flye:
https://github.com/mikolmogorov/Flye/blob/flye/docs/USAGE.md

not-used-snippets/ : moved scripts not used or combined.

about log

log level 0: user interaction or main steps, e.g., skipping ...
log level 1: I/O for the main steps
log level 2: details of I/O
log level 3: function-level and command executions

# 2024-10-13

## Sequence length from Fasta

https://www.biostars.org/p/118954/

# 2024-10-02

## Compare POLAP vith PMAT

conda install conda-forge::apptainer
conda install bioconda::canu
conda install bioconda::nextdenovo

https://nextdenovo.readthedocs.io/en/latest/OPTION.html

# 2024-05-01 (roughly)

## Semi-auto seed contig selection

conda install networkx
conda install gftools
conda install pandas
conda install r::r-mixr

# FIXME

╰─ grep \_polap_var_ga_annotation_all polap-variables-\*.sh ─╯
polap-variables-mtcontig.sh:local \_polap_var_ga_annotation_all="${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt"
polap-variables-mtdna.sh:local _polap_var_ga_annotation_all="${\_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt"

╰─ grep \_polap_var_ga_contigger_edges_gfa polap-variables-\*.sh ─╯
polap-variables-ga.sh:local \_polap_var_ga_contigger_edges_gfa="${_polap_var_ga_contigger}/graph_final.gfa"
polap-variables-mtcontig.sh:local _polap_var_ga_contigger_edges_gfa="${\_polap_var_ga_contigger}/graph_final.gfa"

╰─ grep \_polap_var_wga_annotation polap-variables-\*.sh ─╯
polap-variables-ga.sh:local \_polap_var_wga_annotation="${_polap_var_wga}/assembly_info_organelle_annotation_count-all.txt"
polap-variables-wga.sh:local _polap_var_wga_annotation="${\_polap_var_wga}/assembly_info_organelle_annotation_count-all.txt"

╰─ grep \_polap_var_wga_contigger polap-variables-\*.sh ─╯
polap-variables-ga.sh:local \_polap_var_wga_contigger="${_polap_var_wga}/30-contigger"
polap-variables-ga.sh:local _polap_var_wga_contigger_edges_gfa="${\_polap_var_wga_contigger}/graph_final.gfa"
polap-variables-wga.sh:local \_polap_var_wga_contigger="${_polap_var_wga}/30-contigger"
polap-variables-wga.sh:local _polap_var_wga_contigger_edges_gfa="${\_polap_var_wga_contigger}/graph_final.gfa"
polap-variables-wga.sh:local \_polap_var_wga_contigger_contigs_stats="${_polap_var_wga_contigger}/contigs_stats.txt"
polap-variables-wga.sh:local _polap_var_wga_contigger_contigs_fasta="${\_polap_var_wga_contigger}/contigs.fasta"
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

# # base folders
# local _polap_var_manual_depth_range="${_polap_var_ga}/1-manual.depth.range.txt"
# local _polap_var_manual_copy_range="${_polap_var_ga}/1-manual.copy.range.txt"

# local _polap_var_1_custom_depth_range="${_polap_var_ga}/1-custom.depth.range.txt"
# local _polap_var_2_custom_depth_range="${_polap_var_ga}/2-custom.depth.range.txt"
# local _polap_var_2_depth_range_by_cdf_copy_number="${_polap_var_ga}/2-depth.range.by.cdf.copy.number.txt"
# local _polap_var_3_depth_range_by_mixture="${_polap_var_ga}/3-depth.range.by.mixture.txt"
# local _polap_var_3_mixfit="${_polap_var_ga}/3-mixfit.txt"


# # mtcontigs
# local _polap_var_mtcontigs_1_custom_depth_range="${_polap_var_mtcontigs}/1-custom.depth.range.txt"
# local _polap_var_mtcontigs_2_custom_depth_range="${_polap_var_mtcontigs}/2-custom.depth.range.txt"
# local _polap_var_mtcontigs_2_depth_range_by_cdf_copy_number="${_polap_var_mtcontigs}/2-depth.range.by.cdf.copy.number.txt"
# local _polap_var_mtcontigs_3_depth_range_by_mixture="${_polap_var_mtcontigs}/3-depth.range.by.mixture.txt"
# local _polap_var_mtcontigs_3_mixfit="${_polap_var_mtcontigs}/3-mixfit.txt"

# # output
# local _polap_var_mtcontig_table="${_polap_var_mtcontigs}/1-mtcontig.table.tsv"

# # run-polap-r-select-contigs-by-1-preselection.R
# local _polap_var_preselection_by_gene_density="${_polap_var_mtcontigs}/1-preselection.by.gene.density.txt"
# local _polap_var_preselection_by_depth_mixture="${_polap_var_mtcontigs}/1-preselection.by.depth.mixture.txt"
# local _polap_var_depth_range_by_gene_density="${_polap_var_mtcontigs}/1-depth.range.by.gene.density.txt"
# local _polap_var_depth_range_by_depth_mixture="${_polap_var_mtcontigs}/1-depth.range.by.depth.mixture.txt"
# local _polap_var_mixfit="${_polap_var_mtcontigs}/1-mixfit.txt"

# # run-polap-r-determine-depth-range.R
# local _polap_var_depth_range_by_cdf_copy_number="${_polap_var_mtcontigs}/2-depth.range.by.cdf.copy.number.txt"

# # run-polap-select-contigs-by-1-annotation.R
# local _polap_var_mtcontig_depth_range="${_polap_var_mtcontigs}/1-depth.range.txt"

# local _polap_var_mtcontig_base="${_polap_var_mtcontigs}/1-mtcontig"
# local _polap_var_mtcontig_annotated="${_polap_var_mtcontigs}/1-mtcontig.annotated.txt"
# local _polap_var_mtcontig_annotated_stats="${_polap_var_mtcontigs}/1-mtcontig.annotated.stats.txt"

# # Step 3: filtering GFA using the depth range from step 2
# # run-polap-select-contigs-by-3-subset-by-depth.R
# local _polap_var_mtcontigs_gfa_all="${_polap_var_mtcontigs}/3-gfa.all.gfa"
# local _polap_var_mtcontigs_gfa_seq_part="${_polap_var_mtcontigs}/3-gfa.seq.part.tsv"
# local _polap_var_mtcontigs_gfa_seq_filtered="${_polap_var_mtcontigs}/3-gfa.seq.filtered.txt"

# # Step 4: preparing the graph for finding connected components
# local _polap_var_mtcontigs_gfa_seq_filtered_edge="${_polap_var_mtcontigs}/2-gfa.seq.filtered.edge.txt"
# local _polap_var_mtcontigs_gfa_depthfiltered_gfa="${_polap_var_mtcontigs}/2-gfa.filtered.gfa"

# local _polap_var_mtcontigs_links_tsv="${_polap_var_mtcontigs_links}.tsv"
# local _polap_var_mtcontigs_links_number="${_polap_var_mtcontigs_links}.number.txt"
# local _polap_var_mtcontigs_links_order="${_polap_var_mtcontigs_links}.order.txt"
# local _polap_var_mtcontigs_links_contig="${_polap_var_mtcontigs_links}.contig.txt"
# local _polap_var_mtcontigs_links_contig_na="${_polap_var_mtcontigs_links}.contig.na.txt"

# # Step 5: Find connected components using Python script
# local _polap_var_mtcontigs_links_seed="${_polap_var_mtcontigs}/5-gfa.links.seed.txt"

# # Step 6: converting the depth-filtered contigs in edge with numbers
# local _polap_var_mtcontigs_links_mtcontig="${_polap_var_mtcontigs}/6-gfa.links.mtcontig.txt"
