#/usr/bin/bash

ODIR=ODIR
INUM=INUM
JNUM=JNUM
KNUM=KNUM

script_dir=src

common() {
	source "$script_dir/polap-variables-common.sh"
	declare -p \
		_polap_var_output \
		_polap_var_base \
		_polap_var_project \
		_polap_var_wga \
		_polap_var_ga \
		_polap_var_contigger \
		_polap_var_ga_contigger \
		_polap_var_ann \
		MTCONTIGNAME \
		_polap_var_mtcontigname \
		_polap_var_oga \
		_polap_var_oga_contigger \
		_polap_var_oga \
		_polap_var_oga_seeds \
		_polap_var_mtdna \
		_polap_var_compare \
		_polap_var_base_fq_stats \
		_polap_var_base_genome_size \
		_polap_var_base_jellyfish_out \
		_polap_var_base_jellyfish_out_histo \
		_polap_var_base_l_fq_gz \
		_polap_var_base_long_total_length \
		_polap_var_base_msbwt \
		_polap_var_base_msbwt_tar_gz \
		_polap_var_base_nk_fq_gz \
		_polap_var_base_nk_fq_stats \
		_polap_var_base_lk_fq_gz \
		_polap_var_base_lk_fq_stats \
		_polap_var_bioproject_runinfo_all \
		_polap_var_bioproject_txt \
		_polap_var_bioproject_blastn1 \
		_polap_var_bioproject_blastn2 \
		_polap_var_bioproject_blastn3 \
		_polap_var_bioproject_blastn3_length \
		_polap_var_bioproject_mtdna_fasta1 \
		_polap_var_bioproject_mtdna_fasta1_stats \
		_polap_var_bioproject_mtdna_fasta2 \
		_polap_var_bioproject_mtdna_fasta2_accession \
		_polap_var_bioproject_passed \
		_polap_var_bioproject_runinfo \
		_polap_var_bioproject_species \
		_polap_var_bioproject_sra_long_read \
		_polap_var_bioproject_sra_per_species \
		_polap_var_bioproject_sra_short_read \
		_polap_var_bioproject_taxon_id \
		_polap_var_bioproject_taxonomy \
		_polap_var_wga_annotation \
		_polap_var_wga_contigger \
		_polap_var_wga_contigger_contigs_fasta \
		_polap_var_wga_contigger_contigs_stats \
		_polap_var_wga_contigger_gfa \
		_polap_var_annotation_table \
		_polap_var_ga_annotation_all \
		_polap_var_ga_annotation \
		_polap_var_ga_annotation_depth_table \
		_polap_var_ga_annotation_table \
		_polap_var_ga_annotation_table_contig \
		_polap_var_ga_annotation_all_backup \
		_polap_var_ga_annotation_cdf_table \
		_polap_var_ga_annotation_depth_table_seed \
		_polap_var_ga_annotation_depth_table_seed_target \
		_polap_var_ga_pt_annotation_depth_table \
		_polap_var_ga_pt_annotation_depth_table_seed_target \
		_polap_var_ga_mtcontigs \
		_polap_var_manual_copy_range \
		_polap_var_manual_depth_range \
		_polap_var_assembly_graph_final_gfa \
		_polap_var_contigger_contigs_fasta \
		_polap_var_contigger_contigs_stats \
		_polap_var_contigger_edges_fasta \
		_polap_var_contigger_edges_gfa \
		_polap_var_contigger_edges_stats \
		_polap_var_contigger_gfa \
		_polap_var_ga_contigger_contigs_fasta \
		_polap_var_ga_contigger_contigs_stats \
		_polap_var_ga_contigger_edges_fasta \
		_polap_var_ga_contigger_edges_gfa \
		_polap_var_ga_contigger_edges_stats \
		_polap_var_ga_contigger_gfa \
		_polap_var_ga_gfa_all \
		_polap_var_ga_gfa_seq_part \
		_polap_var_ann_CONTIGDB \
		_polap_var_ann_CONTIGFILE \
		_polap_var_ann_CONTIGNAME \
		_polap_var_ann_MTAABED \
		_polap_var_ann_MTAABLASTBED \
		_polap_var_ann_MTAABLAST \
		_polap_var_ann_MTGENECOUNT \
		_polap_var_ann_MTGENECOUNT_TMP \
		_polap_var_ann_PTAABED \
		_polap_var_ann_PTAABLASTBED \
		_polap_var_ann_PTAABLAST \
		_polap_var_ann_PTGENECOUNT \
		_polap_var_ann_PTGENECOUNT_TMP \
		_polap_var_oga_assembly_graph_gfa \
		_polap_var_oga_contigger_edges_fasta \
		_polap_var_oga_contigger_edges_gfa \
		_polap_var_compare_mtdna1 \
		_polap_var_compare_mtdna2 \
		_polap_var_compare_mtdna3 \
		_polap_var_compare_mtdna_compare \
		_polap_var_compare_oga_blastn1 \
		_polap_var_compare_oga_blastn2 \
		_polap_var_compare_oga_blastn3_length \
		_polap_var_compare_oga_blastn3 \
		_polap_var_oga_contig \
		_polap_var_oga_reads \
		_polap_var_oga_seeds \
		_polap_var_oga_sample \
		_polap_var_oga_flye \
		_polap_var_oga_summary \
		_polap_var_oga_plot \
		_polap_var_mtdna_1_gfa_all \
		_polap_var_mtdna_gfa_links \
		_polap_var_mtdna_gfa_links_edges \
		_polap_var_mtdna_gfa_links_circular_path \
		_polap_var_mtdna_circular_path \
		_polap_var_mtdna_edge_fasta \
		_polap_var_mtdna_fasta \
		_polap_var_mt_fasta \
		_polap_var_mt_edges
}

mtcontig() {
	source "$script_dir/polap-variables-common.sh"
	source "$script_dir/polap-variables-mtcontigs.sh"
	declare -p _polap_var_ga_mtcontigs \
		_polap_var_mtcontigs \
		_polap_var_mtcontigs_depth_range_preselection \
		_polap_var_mtcontigs_depth_range_graphfilter \
		_polap_var_mtcontigs_1_custom_depth_range \
		_polap_var_mtcontigs_2_custom_depth_range \
		_polap_var_mtcontigs_preselection \
		_polap_var_mtcontigs_gfa_all \
		_polap_var_mtcontigs_gfa_seq_filtered \
		_polap_var_mtcontigs_gfa_seq_part \
		_polap_var_mtcontigs_gfa_seq_filtered_edge \
		_polap_var_mtcontigs_gfa_depthfiltered_gfa \
		_polap_var_mtcontigs_links \
		_polap_var_mtcontigs_links_tsv \
		_polap_var_mtcontigs_links_contig_na \
		_polap_var_mtcontigs_links_contig \
		_polap_var_mtcontigs_links_number \
		_polap_var_mtcontigs_links_order \
		_polap_var_mtcontigs_links_seed \
		_polap_var_mtcontigs_links_mtcontig \
		_polap_var_mtcontigs_7mtcontigname \
		_polap_var_mtcontigs_8mtcontigname \
		_polap_var_mtcontig_table \
		_polap_var_mtcontigs_annotation_table_seed \
		_polap_var_mtcontigs_2_depth_range_by_cdf_copy_number
}

echo "-----------------------------------------------------------------------"
echo "common:"
common
echo "-----------------------------------------------------------------------"
echo "mtcontig:"
mtcontig
