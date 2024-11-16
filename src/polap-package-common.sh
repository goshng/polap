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
# _arg_outdir
# _arg_inum
# _arg_jnum

# common variables used across multiple scripts
local _ppack_var_output="${_arg_archive}"
local _ppack_var_base="${_ppack_var_output}"
local _ppack_var_project="${_ppack_var_output}/00-bioproject"
local _ppack_var_wga="${_ppack_var_output}/0"
local _ppack_var_ga="${_ppack_var_output}/${_arg_inum}"
local _ppack_var_contigger="${_ppack_var_ga}/30-contigger"
# _polap_var_ga_contigger better than _polap_var_ga_contigger
local _ppack_var_ga_contigger="${_ppack_var_ga}/30-contigger"
local _ppack_var_ann="${_ppack_var_ga}/50-annotation"
local _polap_var_mtcontigname="${_ppack_var_ga}/mt.contig.name-${_arg_jnum}"
local _ppack_var_mtcontigname="${_ppack_var_ga}/mt.contig.name-${_arg_jnum}"
local _ppack_var_oga="${_ppack_var_output}/${_arg_jnum}"
local _ppack_var_oga_contigger="${_ppack_var_oga}/30-contigger"
# local _ppack_var_oga="${_ppack_var_ga}"
local _ppack_var_oga_seeds="${_ppack_var_oga}/seeds"
local _ppack_var_mtdna="${_ppack_var_ga}/52-mtdna"
local _ppack_var_compare="${_ppack_var_ga}/53-compare"

# base
local _ppack_var_base_fq_stats="${_ppack_var_base}/fq.stats"
local _ppack_var_base_genome_size="${_ppack_var_base}/short_expected_genome_size.txt"
local _ppack_var_base_jellyfish_out="${_ppack_var_base}/jellyfish_out"
local _ppack_var_base_jellyfish_out_histo="${_ppack_var_base}/jellyfish_out.histo"
local _ppack_var_base_l_fq_gz="${_ppack_var_base}/l.fq.gz"
local _ppack_var_base_long_total_length="${_ppack_var_base}/long_total_length.txt"
local _ppack_var_base_msbwt="${_ppack_var_base}/msbwt/comp_msbwt.npy"
local _ppack_var_base_msbwt_tar_gz="${_ppack_var_base}/msbwt.tar.gz"
local _ppack_var_base_nk_fq_gz="${_ppack_var_base}/nk.fq.gz"
local _ppack_var_base_nk_fq_stats="${_ppack_var_base}/nk.fq.stats"
local _ppack_var_base_lk_fq_gz="${_ppack_var_base}/lk.fq.gz"
local _ppack_var_base_lk_fq_stats="${_ppack_var_base}/lk.fq.stats"
local _ppack_var_bioproject_runinfo_all="${_ppack_var_base}/bioproject.runinfo"
local _ppack_var_bioproject_txt="${_ppack_var_base}/bioproject.txt"

# project
local _ppack_var_bioproject_blastn1="${_ppack_var_project}/3-blastn1.txt"
local _ppack_var_bioproject_blastn2="${_ppack_var_project}/3-blastn2.txt"
local _ppack_var_bioproject_blastn3="${_ppack_var_project}/3-blastn3.txt"
local _ppack_var_bioproject_blastn3_length="${_ppack_var_project}/3-blastn3.length.txt"
local _ppack_var_bioproject_mtdna_fasta1="${_ppack_var_project}/1-mtdna.fasta"
local _ppack_var_bioproject_mtdna_fasta1_stats="${_ppack_var_project}/1-mtdna.fasta.stats"
local _ppack_var_bioproject_mtdna_fasta2="${_ppack_var_project}/2-mtdna.fasta"
local _ppack_var_bioproject_mtdna_fasta2_accession="${_ppack_var_project}/2-mtdna.accession"
local _ppack_var_bioproject_passed="${_ppack_var_project}/1-passed.txt"
local _ppack_var_bioproject_runinfo="${_ppack_var_project}/1-runinfo.tsv"
local _ppack_var_bioproject_species="${_ppack_var_project}/1-species.txt"
local _ppack_var_bioproject_sra_long_read="${_ppack_var_project}/1-sra-long-read.tsv"
local _ppack_var_bioproject_sra_per_species="${_ppack_var_project}/1-runinfo.per.species.tsv"
local _ppack_var_bioproject_sra_short_read="${_ppack_var_project}/1-sra-short-read.tsv"
local _ppack_var_bioproject_taxon_id="${_ppack_var_project}/1-taxon-id.txt"
local _ppack_var_bioproject_taxonomy="${_ppack_var_project}/1-taxonomy.txt"

# wga
local _ppack_var_wga_annotation="${_ppack_var_wga}/assembly_info_organelle_annotation_count-all.txt"
local _ppack_var_wga_contigger="${_ppack_var_wga}/30-contigger"
local _ppack_var_wga_contigger_contigs_fasta="${_ppack_var_wga_contigger}/contigs.fasta"
local _ppack_var_wga_contigger_contigs_stats="${_ppack_var_wga_contigger}/contigs_stats.txt"
local _ppack_var_wga_contigger_gfa="${_ppack_var_wga_contigger}/graph_final.gfa"

# ga
local _ppack_var_annotation_table="${_ppack_var_ga}/assembly_info_organelle_annotation_count-all.txt"
local _ppack_var_ga_annotation_all="${_ppack_var_ga}/assembly_info_organelle_annotation_count-all.txt"
local _ppack_var_ga_annotation="${_ppack_var_ga}/assembly_info_organelle_annotation_count.txt"
local _ppack_var_ga_annotation_depth_table="${_ppack_var_ga}/contig-annotation-depth-table.txt"
local _ppack_var_ga_annotation_table="${_ppack_var_ga}/contig-annotation-table.txt"
# for v0.2.6; the contig version not edge one
local _ppack_var_ga_annotation_table_contig="${_ppack_var_ga}/contig-annotation-table-contig.txt"
local _ppack_var_ga_annotation_all_backup="${_ppack_var_ga}/assembly_info_organelle_annotation_count-all.txt.backup"
local _ppack_var_ga_annotation_cdf_table="${_ppack_var_ga}/contig-annotation-cdf-table.txt"
local _ppack_var_ga_annotation_depth_table_seed="${_ppack_var_ga}/contig-annotation-depth-table-seed.txt"
local _ppack_var_ga_annotation_depth_table_seed_target="${_ppack_var_ga}/contig-annotation-depth-table-seed-${_arg_jnum}.txt"
# for PT annotation table
local _ppack_var_ga_pt_annotation_depth_table="${_ppack_var_ga}/pt-contig-annotation-depth-table.txt"
local _ppack_var_ga_pt_annotation_depth_table_seed_target="${_ppack_var_ga}/pt-contig-annotation-depth-table-seed-${_arg_jnum}.txt"
# local _ppack_var_ga_pt_annotation_table="${_ppack_var_ga}/pt-contig-annotation-table.txt"

# not used any more: used to be a folder for _polap_var_ga_gfa_all _polap_var_ga_gfa_seq_part
# FIXME: delete these three lines
local _ppack_var_ga_mtcontigs="${_ppack_var_ga}/mtcontigs"
local _ppack_var_manual_copy_range="${_ppack_var_ga}/1-manual.copy.range.txt"
local _ppack_var_manual_depth_range="${_ppack_var_ga}/1-manual.depth.range.txt"

# _polap_var_ga_contigger
local _ppack_var_assembly_graph_final_gfa="${_ppack_var_contigger}/graph_final.gfa"
local _ppack_var_contigger_contigs_fasta="${_ppack_var_contigger}/contigs.fasta"
local _ppack_var_contigger_contigs_stats="${_ppack_var_contigger}/contigs_stats.txt"
local _ppack_var_contigger_edges_fasta="${_ppack_var_contigger}/graph_final.fasta"
local _ppack_var_contigger_edges_gfa="${_ppack_var_contigger}/graph_final.gfa"
local _ppack_var_contigger_edges_stats="${_ppack_var_contigger}/edges_stats.txt"
local _ppack_var_contigger_gfa="${_ppack_var_contigger}/graph_final.gfa"
# ga version of the above
local _ppack_var_ga_contigger_contigs_fasta="${_ppack_var_ga_contigger}/contigs.fasta"
local _ppack_var_ga_contigger_contigs_stats="${_ppack_var_ga_contigger}/contigs_stats.txt"
local _ppack_var_ga_contigger_edges_fasta="${_ppack_var_ga_contigger}/graph_final.fasta"
local _ppack_var_ga_contigger_edges_gfa="${_ppack_var_ga_contigger}/graph_final.gfa"
local _ppack_var_ga_contigger_edges_stats="${_ppack_var_ga_contigger}/edges_stats.txt"
local _ppack_var_ga_contigger_gfa="${_ppack_var_ga_contigger}/graph_final.gfa"
#
local _ppack_var_ga_gfa_all="${_ppack_var_contigger}/3-gfa.all.gfa"
local _ppack_var_ga_gfa_seq_part="${_ppack_var_contigger}/3-gfa.seq.part.tsv"

# ann
local _ppack_var_ann_CONTIGDB="${_ppack_var_ann}/contig"
local _ppack_var_ann_CONTIGFILE="${_ppack_var_ann}/contig.fasta"
local _ppack_var_ann_CONTIGNAME="${_ppack_var_ann}/contig.name"
local _ppack_var_ann_MTAABED="${_ppack_var_ann}/mtaa.bed"
local _ppack_var_ann_MTAABLASTBED="${_ppack_var_ann}/mtaa.blast.bed"
local _ppack_var_ann_MTAABLAST="${_ppack_var_ann}/mtaa.blast"
local _ppack_var_ann_MTGENECOUNT="${_ppack_var_ann}/mt.gene.count"
local _ppack_var_ann_MTGENECOUNT_TMP="${_ppack_var_ann}/mt.gene.count.tmp"
local _ppack_var_ann_PTAABED="${_ppack_var_ann}/ptaa.bed"
local _ppack_var_ann_PTAABLASTBED="${_ppack_var_ann}/ptaa.blast.bed"
local _ppack_var_ann_PTAABLAST="${_ppack_var_ann}/ptaa.blast"
local _ppack_var_ann_PTGENECOUNT="${_ppack_var_ann}/pt.gene.count"
local _ppack_var_ann_PTGENECOUNT_TMP="${_ppack_var_ann}/pt.gene.count.tmp"

# common: mtcontigs
# moved to polap-variables-mtcontigs.sh

# oga
local _ppack_var_oga_assembly_graph_gfa="${_ppack_var_oga}/assembly_graph.gfa"
local _ppack_var_oga_contigger_edges_fasta="${_ppack_var_oga_contigger}/graph_final.fasta"
local _ppack_var_oga_contigger_edges_gfa="${_ppack_var_oga_contigger}/graph_final.gfa"

# in oga: _polap_var_compare
local _ppack_var_compare_mtdna1="${_ppack_var_compare}/mt.1.fa"
local _ppack_var_compare_mtdna2="${_ppack_var_compare}/mt.2.fa"
local _ppack_var_compare_mtdna3="${_ppack_var_compare}/mt.3.fa"
local _ppack_var_compare_mtdna_compare="${_ppack_var_compare}/mt.compare.txt"
local _ppack_var_compare_oga_blastn1="${_ppack_var_compare}/3-blastn1.txt"
local _ppack_var_compare_oga_blastn2="${_ppack_var_compare}/3-blastn2.txt"
local _ppack_var_compare_oga_blastn3_length="${_ppack_var_compare}/3-blastn3.length.txt"
local _ppack_var_compare_oga_blastn3="${_ppack_var_compare}/3-blastn3.txt"

# seeds
local _ppack_var_oga_contig="${_ppack_var_oga}/01-contig"
local _ppack_var_oga_reads="${_ppack_var_oga}/02-reads"
local _ppack_var_oga_seeds="${_ppack_var_oga}/03-seeds"
local _ppack_var_oga_sample="${_ppack_var_oga}/04-subsample"
local _ppack_var_oga_flye="${_ppack_var_oga}/05-flye"
local _ppack_var_oga_summary="${_ppack_var_oga}/06-summary"
local _ppack_var_oga_plot="${_ppack_var_oga}/07-plot"

# mtdna
local _ppack_var_mtdna_1_gfa_all="${_ppack_var_mtdna}/1-gfa.all.gfa"
local _ppack_var_mtdna_gfa_links="${_ppack_var_mtdna}/1-gfa.links.tsv"
local _ppack_var_mtdna_gfa_links_edges="${_ppack_var_mtdna}/1-gfa.links.edges.txt"
local _ppack_var_mtdna_gfa_links_circular_path="${_ppack_var_mtdna}/2-gfa.links.circular.path.txt"
local _ppack_var_mtdna_circular_path="${_ppack_var_mtdna}/3-circular.path.txt"
local _ppack_var_mtdna_edge_fasta="${_ppack_var_mtdna}/5-edge.fasta"
local _ppack_var_mtdna_fasta="${_ppack_var_mtdna}/4-gfa.fasta"

# mt
local _ppack_var_mt_fasta="${_ppack_var_ga}/mt.0.fasta"
local _ppack_var_mt_edges="${_ppack_var_ga}/mt.0.edges"
