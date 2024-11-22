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
local _polap_var_outdir="${_arg_outdir}"
local _polap_var_project="${_polap_var_outdir}/00-bioproject"
local _polap_var_taxonomy="${_polap_var_outdir}/01-taxonomy"
local _polap_var_wga="${_polap_var_outdir}/0"
local _polap_var_ga="${_polap_var_outdir}/${_arg_inum}"
local _polap_var_ga_contigger="${_polap_var_ga}/30-contigger"
# _polap_var_ga_contigger better than _polap_var_ga_contigger
local _polap_var_ga_contigger="${_polap_var_ga}/30-contigger"
local _polap_var_ann="${_polap_var_ga}/50-annotation"
local _polap_var_mtcontigname="${_polap_var_ga}/mt.contig.name-${_arg_jnum}"
local _polap_var_oga="${_polap_var_outdir}/${_arg_jnum}"
local _polap_var_oga_contigger="${_polap_var_oga}/30-contigger"
# local _polap_var_oga="${_polap_var_ga}"
local _polap_var_oga_seeds="${_polap_var_oga}/seeds"
local _polap_var_mtdna="${_polap_var_oga}/52-mtdna"
local _polap_var_compare="${_polap_var_oga}/53-compare"
local _polap_var_geseq="${_polap_var_outdir}/01-geseq"

# base
local _polap_var_outdir_s1_fq_stats="${_polap_var_outdir}/s1.fq.stats"
local _polap_var_outdir_s2_fq_stats="${_polap_var_outdir}/s2.fq.stats"
local _polap_var_outdir_genome_size="${_polap_var_outdir}/short_expected_genome_size.txt"
local _polap_var_outdir_jellyfish_out="${_polap_var_outdir}/jellyfish_out"
local _polap_var_outdir_jellyfish_out_histo="${_polap_var_outdir}/jellyfish_out.histo"
local _polap_var_outdir_l_fq_gz="${_polap_var_outdir}/l.fq.gz"
local _polap_var_outdir_l_fq_stats="${_polap_var_outdir}/l.fq.stats"
local _polap_var_outdir_long_total_length="${_polap_var_outdir}/long_total_length.txt"
local _polap_var_outdir_msbwt_dir="${_polap_var_outdir}/msbwt"
local _polap_var_outdir_msbwt="${_polap_var_outdir}/msbwt/comp_msbwt.npy"
local _polap_var_outdir_msbwt_tar_gz="${_polap_var_outdir}/msbwt.tar.gz"
local _polap_var_outdir_nk_fq_gz="${_polap_var_outdir}/nk.fq.gz"
local _polap_var_outdir_nk_fq_stats="${_polap_var_outdir}/nk.fq.stats"
local _polap_var_outdir_lk_fq_gz="${_polap_var_outdir}/lk.fq.gz"
local _polap_var_outdir_lk_fq_stats="${_polap_var_outdir}/lk.fq.stats"
local _polap_var_project_runinfo_all="${_polap_var_outdir}/bioproject.runinfo"
local _polap_var_project_txt="${_polap_var_outdir}/bioproject.txt"

# project
local _polap_var_project_blastn1="${_polap_var_project}/3-blastn1.txt"
local _polap_var_project_blastn2="${_polap_var_project}/3-blastn2.txt"
local _polap_var_project_blastn3="${_polap_var_project}/3-blastn3.txt"
local _polap_var_project_blastn3_length="${_polap_var_project}/3-blastn3.length.txt"
local _polap_var_project_mtdna_fasta1="${_polap_var_project}/1-mtdna.fasta"
local _polap_var_project_mtdna_fasta1_stats="${_polap_var_project}/1-mtdna.fasta.stats"
local _polap_var_project_mtdna_fasta2="${_polap_var_project}/2-mtdna.fasta"
local _polap_var_project_mtdna_fasta2_accession="${_polap_var_project}/2-mtdna.accession"
local _polap_var_project_passed="${_polap_var_project}/1-passed.txt"
local _polap_var_project_runinfo="${_polap_var_project}/1-runinfo.tsv"
local _polap_var_project_species="${_polap_var_project}/1-species.txt"
local _polap_var_project_sra_long_read="${_polap_var_project}/1-sra-long-read.tsv"
local _polap_var_project_sra_per_species="${_polap_var_project}/1-runinfo.per.species.tsv"
local _polap_var_project_sra_short_read="${_polap_var_project}/1-sra-short-read.tsv"
local _polap_var_project_taxon_id="${_polap_var_project}/1-taxon-id.txt"
local _polap_var_project_taxonomy="${_polap_var_project}/1-taxonomy.txt"

# wga
local _polap_var_wga_annotation_all="${_polap_var_wga}/assembly_info_organelle_annotation_count-all.txt"
local _polap_var_wga_contigger="${_polap_var_wga}/30-contigger"
local _polap_var_wga_contigger_contigs_fasta="${_polap_var_wga_contigger}/contigs.fasta"
local _polap_var_wga_contigger_contigs_stats="${_polap_var_wga_contigger}/contigs_stats.txt"
local _polap_var_wga_contigger_edges_gfa="${_polap_var_wga_contigger}/graph_final.gfa"
local _polap_var_wga_contigger_edges_fasta="${_polap_var_wga_contigger}/graph_final.fasta"
local _polap_var_wga_contigger_edges_stats="${_polap_var_wga_contigger}/edges_stats.txt"

# ga
local _polap_var_ga_annotation_all="${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt"
local _polap_var_ga_annotation="${_polap_var_ga}/assembly_info_organelle_annotation_count.txt"
local _polap_var_ga_annotation_depth_table="${_polap_var_ga}/contig-annotation-depth-table.txt"
local _polap_var_ga_annotation_table="${_polap_var_ga}/contig-annotation-table.txt"
# for v0.2.6; the contig version not edge one
local _polap_var_ga_annotation_table_contig="${_polap_var_ga}/contig-annotation-table-contig.txt"
local _polap_var_ga_annotation_all_backup="${_polap_var_ga}/assembly_info_organelle_annotation_count-all.txt.backup"
local _polap_var_ga_annotation_cdf_table="${_polap_var_ga}/contig-annotation-cdf-table.txt"
# delete table_seed not table_seed_target
local _polap_var_ga_annotation_depth_table_seed="${_polap_var_ga}/contig-annotation-depth-table-seed.txt"
local _polap_var_ga_annotation_depth_table_seed_target="${_polap_var_ga}/contig-annotation-depth-table-seed-${_arg_jnum}.txt"
# for PT annotation table
local _polap_var_ga_pt_annotation_depth_table="${_polap_var_ga}/pt-contig-annotation-depth-table.txt"
local _polap_var_ga_pt_annotation_depth_table_seed_target="${_polap_var_ga}/pt-contig-annotation-depth-table-seed-${_arg_jnum}.txt"

# _polap_var_ga_contigger
local _polap_var_ga_contigger_edges_gfa="${_polap_var_ga_contigger}/graph_final.gfa"
local _polap_var_ga_contigger_contigs_fasta="${_polap_var_ga_contigger}/contigs.fasta"
local _polap_var_ga_contigger_contigs_stats="${_polap_var_ga_contigger}/contigs_stats.txt"
local _polap_var_ga_contigger_edges_fasta="${_polap_var_ga_contigger}/graph_final.fasta"
local _polap_var_ga_contigger_edges_stats="${_polap_var_ga_contigger}/edges_stats.txt"
#
local _polap_var_ga_gfa_all="${_polap_var_ga_contigger}/3-gfa.all.gfa"
local _polap_var_ga_gfa_seq_part="${_polap_var_ga_contigger}/3-gfa.seq.part.tsv"

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
# moved to polap-variables-mtcontigs.sh

# ga
local _polap_var_ga_assembly_graph_gfa="${_polap_var_ga}/assembly_graph.gfa"
local _polap_var_ga_assembly_graph_edges_fasta="${_polap_var_ga}/edges_graph.fasta"
local _polap_var_ga_assembly_graph_edges_stats="${_polap_var_ga}/edges_stats.txt"

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
local _polap_var_oga_subsample="${_polap_var_oga}/04-subsample"
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
local _polap_var_mtdna_names="${_polap_var_mtdna}/1-mtdna.names"
local _polap_var_mtdna_sub_gfa="${_polap_var_mtdna}/9-gfa"

# mt
local _polap_var_mt_fasta="${_polap_var_ga}/mt.0.fasta"
local _polap_var_mt_edges="${_polap_var_ga}/mt.0.edges"
local _polap_var_oga_mt_fasta="${_polap_var_oga}/mt.0.fasta"
local _polap_var_oga_mt_edges="${_polap_var_oga}/mt.0.edges"

local _polap_var_mt1fa="${_polap_var_outdir}/mt.1.fa"
local _polap_var_geseq_gff="${_polap_var_geseq}/mt.1.gff"
local _polap_var_geseq_faa="${_polap_var_geseq}/mt.1.faa"
local _polap_var_geseq_fna="${_polap_var_geseq}/mt.1.fna"
local _polap_var_geseq_filtered_faa="${_polap_var_geseq}/mt.1.filtered.faa"
local _polap_var_geseq_formated_faa="${_polap_var_geseq}/mt.1.formated.faa"

local _polap_var_orthofinder="${_polap_var_outdir}/orthofinder"

# ncbi
local _polap_var_ncbi="${_polap_var_outdir}/ncbi"
local _polap_var_ncbi_accessions="${_polap_var_ncbi}/accession.txt"
local _polap_var_ncbi_sequence_tsv="${_polap_var_ncbi}/sequence.tsv"
local _polap_var_ncbi_sequence_taxon_tsv="${_polap_var_ncbi}/sequence.taxon.tsv"
local _polap_var_ncbi_sequence_taxon_map="${_polap_var_ncbi}/sequence.taxon.map"
local _polap_var_ncbi_sampled_accession="${_polap_var_ncbi}/sampled.accesssion.txt"

# tree
local _polap_var_phylogenies="${_polap_var_outdir}/phylogenies"
local _polap_var_phylogenies_tree="${_polap_var_phylogenies}/s.treefile"
local _polap_var_phylogenies_tree_taxa="${_polap_var_phylogenies}/taxa.treefile"

# output
local _polap_var_output_wga_gfa="1-wga.gfa"
local _polap_var_output_oga_gfa="2-oga.gfa"
