#/usr/bin/bash

ODIR=ODIR
INUM=INUM
JNUM=JNUM
script_dir=src
common() {

	source "$script_dir/polap-variables-common.sh"
	# source "$script_dir/polap-variables-main.sh"
	# source "$script_dir/polap-variables-bioproject.sh"
	# source "$script_dir/polap-variables-ga.sh"
	# source "$script_dir/polap-variables-wga.sh"
	# source "$script_dir/polap-variables-oga.sh"
	# source "$script_dir/polap-variables-mtcontig.sh"

	declare -p _polap_var_output _polap_var_base _polap_var_bioproject _polap_var_wga _polap_var_ga _polap_var_contigger _polap_var_ann _polap_var_mtcontigs _polap_var_mtdna _polap_var_compare MTCONTIGNAME _polap_var_mtcontigname _polap_var_links _polap_var_oga _polap_var_oga_contigger _polap_var_seeds
}

base() {

	source "$script_dir/polap-variables-base.sh"
	# source "$script_dir/polap-variables-main.sh"
	# source "$script_dir/polap-variables-bioproject.sh"
	# source "$script_dir/polap-variables-ga.sh"
	# source "$script_dir/polap-variables-wga.sh"
	# source "$script_dir/polap-variables-oga.sh"
	# source "$script_dir/polap-variables-mtcontig.sh"

	declare -p _polap_var_base_jellyfish_out _polap_var_base_jellyfish_out_histo _polap_var_base_genome_size _polap_var_base_long_total_length _polap_var_base_nk_fq_gz _polap_var_base_l_fq_gz _polap_var_base_nk_fq_stats _polap_var_base_fq_stats _polap_var_base_msbwt _polap_var_base_msbwt_tar_gz _polap_var_bioproject_txt _polap_var_bioproject_runinfo_all
}

bioproject() {

	# source "$script_dir/polap-variables-base.sh"
	# source "$script_dir/polap-variables-main.sh"
	source "$script_dir/polap-variables-bioproject.sh"
	# source "$script_dir/polap-variables-ga.sh"
	# source "$script_dir/polap-variables-wga.sh"
	# source "$script_dir/polap-variables-oga.sh"
	# source "$script_dir/polap-variables-mtcontig.sh"

	declare -p _polap_var_bioproject_runinfo _polap_var_bioproject_sra_per_species _polap_var_bioproject_sra_long_read _polap_var_bioproject_sra_short_read _polap_var_bioproject_species _polap_var_bioproject_taxon_id _polap_var_bioproject_taxonomy _polap_var_bioproject_passed _polap_var_bioproject_mtdna_fasta1 _polap_var_bioproject_mtdna_fasta1_stats _polap_var_bioproject_mtdna_fasta2 _polap_var_bioproject_mtdna_fasta2_accession _polap_var_bioproject_blastn1 _polap_var_bioproject_blastn2 _polap_var_bioproject_blastn3 _polap_var_bioproject_blastn3_length
}

ga() {

	# source "$script_dir/polap-variables-base.sh"
	# source "$script_dir/polap-variables-main.sh"
	# source "$script_dir/polap-variables-bioproject.sh"
	source "$script_dir/polap-variables-ga.sh"
	# source "$script_dir/polap-variables-wga.sh"
	# source "$script_dir/polap-variables-oga.sh"
	# source "$script_dir/polap-variables-mtcontig.sh"

	declare -p _polap_var_ga_annotation _polap_var_ga_annotation_all _polap_var_ga_annotation_all_backup _polap_var_ga_annotation_table _polap_var_ga_annotation_depth_table _polap_var_ga_annotation_depth_table_seed _polap_var_ga_annotation_depth_table_seed_target _polap_var_contigger_gfa _polap_var_contigger_edges_stats _polap_var_contigger_edges_fasta _polap_var_contigger_edges_gfa _polap_var_contigger_contigs_stats _polap_var_contigger_contigs_fasta _polap_var_assembly_graph_final_gfa _polap_var_wga_contigger _polap_var_wga_contigger_gfa _polap_var_wga_annotation _polap_var_ga_mtcontigs _polap_var_ga_gfa_all _polap_var_ga_gfa_seq_part _polap_var_ann_CONTIGFILE _polap_var_ann_CONTIGDB _polap_var_ann_CONTIGNAME _polap_var_ann_MTAABLAST _polap_var_ann_MTAABLASTBED _polap_var_ann_MTAABED _polap_var_ann_MTGENECOUNT_TMP _polap_var_ann_MTGENECOUNT _polap_var_ann_PTAABLAST _polap_var_ann_PTAABLASTBED _polap_var_ann_PTAABED _polap_var_ann_PTGENECOUNT_TMP _polap_var_ann_PTGENECOUNT
}

wga() {

	# source "$script_dir/polap-variables-base.sh"
	# source "$script_dir/polap-variables-main.sh"
	# source "$script_dir/polap-variables-bioproject.sh"
	# source "$script_dir/polap-variables-ga.sh"
	source "$script_dir/polap-variables-wga.sh"
	# source "$script_dir/polap-variables-oga.sh"
	# source "$script_dir/polap-variables-mtcontig.sh"

	declare -p _polap_var_wga_annotation _polap_var_wga_contigger _polap_var_wga_contigger_gfa _polap_var_wga_contigger_contigs_stats _polap_var_wga_contigger_contigs_fasta
}

oga() {

	# source "$script_dir/polap-variables-base.sh"
	# source "$script_dir/polap-variables-main.sh"
	# source "$script_dir/polap-variables-bioproject.sh"
	# source "$script_dir/polap-variables-ga.sh"
	# source "$script_dir/polap-variables-wga.sh"
	source "$script_dir/polap-variables-oga.sh"
	# source "$script_dir/polap-variables-mtcontig.sh"

	declare -p _polap_var_oga_contigger_edges_gfa _polap_var_oga_contigger_edges_fasta _polap_var_oga_assembly_graph_gfa _polap_var_mtdna1 _polap_var_mtdna2 _polap_var_mtdna3 _polap_var_mtdna_compare _polap_var_oga_blastn1 _polap_var_oga_blastn2 _polap_var_oga_blastn3 _polap_var_oga_blastn3_length
}

mtcontig() {

	source "$script_dir/polap-variables-mtcontig.sh"

	declare -p _polap_var_manual_depth_range _polap_var_manual_copy_range _polap_var_1_custom_depth_range _polap_var_2_custom_depth_range _polap_var_2_depth_range_by_cdf_copy_number _polap_var_3_depth_range_by_mixture _polap_var_3_mixfit _polap_var_mtcontigs_1_custom_depth_range _polap_var_mtcontigs_2_custom_depth_range _polap_var_mtcontigs_2_depth_range_by_cdf_copy_number _polap_var_mtcontigs_3_depth_range_by_mixture _polap_var_mtcontigs_3_mixfit _polap_var_assembly_graph_final_gfa _polap_var_annotation_table _polap_var_mtcontig_table _polap_var_preselection_by_gene_density _polap_var_preselection_by_depth_mixture _polap_var_depth_range_by_gene_density _polap_var_depth_range_by_depth_mixture _polap_var_mixfit _polap_var_depth_range_by_cdf_copy_number _polap_var_mtcontig_depth_range _polap_var_mtcontig_base _polap_var_mtcontig_annotated _polap_var_mtcontig_annotated_stats _polap_var_gfa_all _polap_var_gfa_seq_part _polap_var_gfa_seq_filtered _polap_var_gfa_seq_filtered_edge _polap_var_gfa_filtered _polap_var_links_tsv _polap_var_links_number _polap_var_links_order _polap_var_links_contig _polap_var_links_contig_na _polap_var_links_seed _polap_var_links_mtcontig
}

echo "-----------------------------------------------------------------------"
echo "common:"
common
echo "-----------------------------------------------------------------------"
echo "base:"
base
echo "-----------------------------------------------------------------------"
echo "bioproject:"
bioproject
echo "-----------------------------------------------------------------------"
echo "ga:"
ga
echo "-----------------------------------------------------------------------"
echo "wga:"
wga
echo "-----------------------------------------------------------------------"
echo "oga:"
oga
echo "-----------------------------------------------------------------------"
echo "mtcontig:"
mtcontig
