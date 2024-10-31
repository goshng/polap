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

# Global option variables
#
# ODIR
# INUM
# JNUM

# common variables used across multiple scripts
local _polap_var_output="${ODIR}"
local _polap_var_base="${_polap_var_output}"
local _polap_var_project="${_polap_var_output}/00-bioproject"
local _polap_var_wga="${_polap_var_output}/0"
local _polap_var_ga="${_polap_var_output}/${INUM}"
local _polap_var_contigger="${_polap_var_ga}/30-contigger"
local _polap_var_ann="${_polap_var_ga}/50-annotation"
local _polap_var_mtcontigs="${_polap_var_ga}/${JNUM}/mtcontigs"
# local _polap_var_mtcontigs="${_polap_var_ga}/60-mtcontigs"
local MTCONTIGNAME="${_polap_var_ga}/mt.contig.name-${JNUM}"
local _polap_var_mtcontigname="${_polap_var_ga}/mt.contig.name-${JNUM}"
local _polap_var_oga="${_polap_var_output}/${JNUM}"
local _polap_var_oga_contigger="${_polap_var_oga}/30-contigger"
# local _polap_var_oga="${_polap_var_ga}"
local _polap_var_oga_seeds="${_polap_var_oga}/seeds"
local _polap_var_links="${_polap_var_mtcontigs}/4-gfa.links"
local _polap_var_mtdna="${_polap_var_ga}/51-mtdna"
local _polap_var_compare="${_polap_var_ga}/52-compare"

# base
local _polap_var_base_fq_stats="${_polap_var_base}/fq.stats"
local _polap_var_base_genome_size="${_polap_var_base}/short_expected_genome_size.txt"
local _polap_var_base_jellyfish_out="${_polap_var_base}/jellyfish_out"
local _polap_var_base_jellyfish_out_histo="${_polap_var_base}/jellyfish_out.histo"
local _polap_var_base_l_fq_gz="${_polap_var_base}/l.fq.gz"
local _polap_var_base_long_total_length="${_polap_var_base}/long_total_length.txt"
local _polap_var_base_msbwt="${_polap_var_base}/msbwt/comp_msbwt.npy"
local _polap_var_base_msbwt_tar_gz="${_polap_var_base}/msbwt.tar.gz"
local _polap_var_base_nk_fq_gz="${_polap_var_base}/nk.fq.gz"
local _polap_var_base_nk_fq_stats="${_polap_var_base}/nk.fq.stats"
local _polap_var_base_lk_fq_gz="${_polap_var_base}/lk.fq.gz"
local _polap_var_base_lk_fq_stats="${_polap_var_base}/lk.fq.stats"
local _polap_var_bioproject_runinfo_all="${_polap_var_base}/bioproject.runinfo"
local _polap_var_bioproject_txt="${_polap_var_base}/bioproject.txt"

# project
local _polap_var_bioproject_blastn1="${_polap_var_project}/3-blastn1.txt"
local _polap_var_bioproject_blastn2="${_polap_var_project}/3-blastn2.txt"
local _polap_var_bioproject_blastn3="${_polap_var_project}/3-blastn3.txt"
local _polap_var_bioproject_blastn3_length="${_polap_var_project}/3-blastn3.length.txt"
local _polap_var_bioproject_mtdna_fasta1="${_polap_var_project}/1-mtdna.fasta"
local _polap_var_bioproject_mtdna_fasta1_stats="${_polap_var_project}/1-mtdna.fasta.stats"
local _polap_var_bioproject_mtdna_fasta2="${_polap_var_project}/2-mtdna.fasta"
local _polap_var_bioproject_mtdna_fasta2_accession="${_polap_var_project}/2-mtdna.accession"
local _polap_var_bioproject_passed="${_polap_var_project}/1-passed.txt"
local _polap_var_bioproject_runinfo="${_polap_var_project}/1-runinfo.tsv"
local _polap_var_bioproject_species="${_polap_var_project}/1-species.txt"
local _polap_var_bioproject_sra_long_read="${_polap_var_project}/1-sra-long-read.tsv"
local _polap_var_bioproject_sra_per_species="${_polap_var_project}/1-runinfo.per.species.tsv"
local _polap_var_bioproject_sra_short_read="${_polap_var_project}/1-sra-short-read.tsv"
local _polap_var_bioproject_taxon_id="${_polap_var_project}/1-taxon-id.txt"
local _polap_var_bioproject_taxonomy="${_polap_var_project}/1-taxonomy.txt"

# wga
local _polap_var_wga_annotation="${_polap_var_wga}/assembly_info_organelle_annotation_count-all.txt"
local _polap_var_wga_contigger="${_polap_var_wga}/30-contigger"
local _polap_var_wga_contigger_contigs_fasta="${_polap_var_wga_contigger}/contigs.fasta"
local _polap_var_wga_contigger_contigs_stats="${_polap_var_wga_contigger}/contigs_stats.txt"
local _polap_var_wga_contigger_gfa="${_polap_var_wga_contigger}/graph_final.gfa"

# ga
local _polap_var_1_custom_depth_range="${_polap_var_ga}/1-custom.depth.range.txt"
local _polap_var_2_custom_depth_range="${_polap_var_ga}/2-custom.depth.range.txt"
local _polap_var_2_depth_range_by_cdf_copy_number="${_polap_var_ga}/2-depth.range.by.cdf.copy.number.txt"
local _polap_var_3_depth_range_by_mixture="${_polap_var_ga}/3-depth.range.by.mixture.txt"
local _polap_var_3_mixfit="${_polap_var_ga}/3-mixfit.txt"
local _polap_var_annotation_table="${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt"
local _polap_var_ga_annotation="${_polap_var_ga}/assembly_info_organelle_annotation_count.txt"
local _polap_var_ga_annotation_all="${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt"
local _polap_var_ga_annotation_all_backup="${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt.backup"
local _polap_var_ga_annotation_depth_table="${_polap_var_ga}/contig-annotation-depth-table.txt"
local _polap_var_ga_annotation_cdf_table="${_polap_var_ga}/contig-annotation-cdf-table.txt"
local _polap_var_ga_annotation_depth_table_seed="${_polap_var_ga}/contig-annotation-depth-table-seed.txt"
local _polap_var_ga_annotation_depth_table_seed_target="${_polap_var_ga}/contig-annotation-depth-table-seed-${JNUM}.txt"
local _polap_var_ga_annotation_table="${_polap_var_ga}/contig-annotation-table.txt"
# not used any more: used to be a folder for _polap_var_ga_gfa_all _polap_var_ga_gfa_seq_part
local _polap_var_ga_mtcontigs="${_polap_var_ga}/mtcontigs"
local _polap_var_manual_copy_range="${_polap_var_ga}/1-manual.copy.range.txt"
local _polap_var_manual_depth_range="${_polap_var_ga}/1-manual.depth.range.txt"
# not used at the moment: could replace MTCONTIGNAME variable
local _polap_var_mtcontigname="${_polap_var_ga}/mt.contig.name-${JNUM}"

# _polap_var_contigger
local _polap_var_assembly_graph_final_gfa="${_polap_var_contigger}/graph_final.gfa"
local _polap_var_contigger_contigs_fasta="${_polap_var_contigger}/contigs.fasta"
local _polap_var_contigger_contigs_stats="${_polap_var_contigger}/contigs_stats.txt"
local _polap_var_contigger_edges_fasta="${_polap_var_contigger}/graph_final.fasta"
local _polap_var_contigger_edges_gfa="${_polap_var_contigger}/graph_final.gfa"
local _polap_var_contigger_edges_stats="${_polap_var_contigger}/edges_stats.txt"
local _polap_var_contigger_gfa="${_polap_var_contigger}/graph_final.gfa"
local _polap_var_ga_gfa_all="${_polap_var_contigger}/3-gfa.all.gfa"
local _polap_var_ga_gfa_seq_part="${_polap_var_contigger}/3-gfa.seq.part.tsv"

# ann
local _polap_var_ann_CONTIGDB="${_polap_var_ann}/contig"
local _polap_var_ann_CONTIGFILE="${_polap_var_ann}/contig.fasta"
local _polap_var_ann_CONTIGNAME="${_polap_var_ann}/contig.name"
local _polap_var_ann_MTAABED="${_polap_var_ann}/mtaa.bed"
local _polap_var_ann_MTAABLASTBED="${_polap_var_ann}/mtaa.blast.bed"
local _polap_var_ann_MTAABLAST="${_polap_var_ann}/mtaa.blast"
local _polap_var_ann_MTGENECOUNT="${_polap_var_ann}/mt.gene.count"
local _polap_var_ann_MTGENECOUNT_TMP="${_polap_var_ann}/mt.gene.count.tmp"
local _polap_var_ann_PTAABED="${_polap_var_ann}/ptaa.bed"
local _polap_var_ann_PTAABLASTBED="${_polap_var_ann}/ptaa.blast.bed"
local _polap_var_ann_PTAABLAST="${_polap_var_ann}/ptaa.blast"
local _polap_var_ann_PTGENECOUNT="${_polap_var_ann}/pt.gene.count"
local _polap_var_ann_PTGENECOUNT_TMP="${_polap_var_ann}/pt.gene.count.tmp"

# common: mtcontigs
local _polap_var_depth_range_by_cdf_copy_number="${_polap_var_mtcontigs}/2-depth.range.by.cdf.copy.number.txt"
local _polap_var_depth_range_by_depth_mixture="${_polap_var_mtcontigs}/1-depth.range.by.depth.mixture.txt"
local _polap_var_depth_range_by_gene_density="${_polap_var_mtcontigs}/1-depth.range.by.gene.density.txt"
local _polap_var_gfa_all="${_polap_var_mtcontigs}/3-gfa.all.gfa"
local _polap_var_gfa_filtered="${_polap_var_mtcontigs}/2-gfa.filtered.gfa"
local _polap_var_gfa_seq_filtered_edge="${_polap_var_mtcontigs}/2-gfa.seq.filtered.edge.txt"
local _polap_var_gfa_seq_filtered="${_polap_var_mtcontigs}/3-gfa.seq.filtered.txt"
local _polap_var_gfa_seq_part="${_polap_var_mtcontigs}/3-gfa.seq.part.tsv"
local _polap_var_links_contig_na="${_polap_var_mtcontigs}/4-gfa.links.contig.na.txt"
local _polap_var_links_contig="${_polap_var_mtcontigs}/4-gfa.links.contig.txt"
local _polap_var_links_mtcontig="${_polap_var_mtcontigs}/6-gfa.links.mtcontig.txt"
local _polap_var_links_number="${_polap_var_mtcontigs}/4-gfa.links.number.txt"
local _polap_var_links="${_polap_var_mtcontigs}/4-gfa.links"
local _polap_var_links_order="${_polap_var_mtcontigs}/4-gfa.links.order.txt"
local _polap_var_links_seed="${_polap_var_mtcontigs}/5-gfa.links.seed.txt"
local _polap_var_links_tsv="${_polap_var_mtcontigs}/4-gfa.links.tsv"
local _polap_var_mixfit="${_polap_var_mtcontigs}/1-mixfit.txt"
local _polap_var_mtcontig_annotated="${_polap_var_mtcontigs}/1-mtcontig.annotated.txt"
local _polap_var_mtcontig_annotated_stats="${_polap_var_mtcontigs}/1-mtcontig.annotated.stats.txt"
local _polap_var_mtcontig_base="${_polap_var_mtcontigs}/1-mtcontig"
local _polap_var_mtcontig_depth_range="${_polap_var_mtcontigs}/1-depth.range.txt"
local _polap_var_mtcontigs_1_custom_depth_range="${_polap_var_mtcontigs}/1-custom.depth.range.txt"
local _polap_var_mtcontigs_2_custom_depth_range="${_polap_var_mtcontigs}/2-custom.depth.range.txt"
local _polap_var_mtcontigs_2_depth_range_by_cdf_copy_number="${_polap_var_mtcontigs}/2-depth.range.by.cdf.copy.number.txt"
local _polap_var_mtcontigs_3_depth_range_by_mixture="${_polap_var_mtcontigs}/3-depth.range.by.mixture.txt"
local _polap_var_mtcontigs_3_mixfit="${_polap_var_mtcontigs}/3-mixfit.txt"
local _polap_var_mtcontig_table="${_polap_var_mtcontigs}/1-mtcontig.table.tsv"
local _polap_var_preselection_by_depth_mixture="${_polap_var_mtcontigs}/1-preselection.by.depth.mixture.txt"
local _polap_var_preselection_by_gene_density="${_polap_var_mtcontigs}/1-preselection.by.gene.density.txt"

# oga
local _polap_var_oga_assembly_graph_gfa="${_polap_var_oga}/assembly_graph.gfa"
local _polap_var_oga_contigger_edges_fasta="${_polap_var_oga_contigger}/graph_final.fasta"
local _polap_var_oga_contigger_edges_gfa="${_polap_var_oga_contigger}/graph_final.gfa"

# in oga: _polap_var_compare
local _polap_var_compare_mtdna1="${_polap_var_compare}/mt.1.fa"
local _polap_var_compare_mtdna2="${_polap_var_compare}/mt.2.fa"
local _polap_var_compare_mtdna3="${_polap_var_compare}/mt.3.fa"
local _polap_var_compare_mtdna_compare="${_polap_var_compare}/mt.compare.txt"
local _polap_var_compare_oga_blastn1="${_polap_var_compare}/3-blastn1.txt"
local _polap_var_compare_oga_blastn2="${_polap_var_compare}/3-blastn2.txt"
local _polap_var_compare_oga_blastn3_length="${_polap_var_compare}/3-blastn3.length.txt"
local _polap_var_compare_oga_blastn3="${_polap_var_compare}/3-blastn3.txt"

# seeds
local _polap_var_oga_contig="${_polap_var_oga}/01-contig"
local _polap_var_oga_reads="${_polap_var_oga}/02-reads"
local _polap_var_oga_seeds="${_polap_var_oga}/03-seeds"
local _polap_var_oga_sample="${_polap_var_oga}/04-subsample"
local _polap_var_oga_flye="${_polap_var_oga}/05-flye"
local _polap_var_oga_summary="${_polap_var_oga}/06-summary"
local _polap_var_oga_plot="${_polap_var_oga}/07-plot"

# mtdna
local _polap_var_mtdna_1_gfa_all="${_polap_var_mtdna}/1-gfa.all.gfa"
local _polap_var_mtdna_gfa_links="${_polap_var_mtdna}/1-gfa.links.tsv"
local _polap_var_mtdna_gfa_links_edges="${_polap_var_mtdna}/1-gfa.links.edges.txt"
local _polap_var_mtdna_gfa_links_circular_path="${_polap_var_mtdna}/2-gfa.links.circular.path.txt"
local _polap_var_mtdna_circular_path="${_polap_var_mtdna}/3-circular.path.txt"
local _polap_var_mtdna_edge_fasta="${_polap_var_mtdna}/5-edge.fasta"
local _polap_var_mtdna_fasta="${_polap_var_mtdna}/4-gfa.fasta"

# mt
local _polap_var_mt_fasta="${_polap_var_ga}/mt.0.fasta"
local _polap_var_mt_edges="${_polap_var_ga}/mt.0.edges"
