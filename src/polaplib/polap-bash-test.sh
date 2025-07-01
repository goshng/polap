#!/usr/bin/env bash

_POLAPLIB_DIR=./src/polaplib

test-bash-extract-three-edges-of-ptdna() {
	bash "${_POLAPLIB_DIR}"/polap-bash-extract-three-edges-of-ptdna.sh input/assembly_graph.gfa output/extract-three-edges-of-ptdna.txt
	diff output/extract-three-edges-of-ptdna.txt output_expected/extract-three-edges-of-ptdna.txt
}

test-r-depth-distribution() {
	# activate polap conda environment
	Rscript "${_POLAPLIB_DIR}"/run-polap-r-depth-distribution.R \
		-t input/assembly_info_organelle_annotation_count-all.txt \
		-o output/assembly_info_organelle_annotation_count-all.txt.pdf
	echo open output/assembly_info_organelle_annotation_count-all.txt.pdf
}

test-bash-create-depth-file() {
	local _polap_var_mtcontigs=output
	local _polap_var_mtcontigs_depth_range_graphfilter="${_polap_var_mtcontigs}/1-depth.range.txt"
	bash "${_POLAPLIB_DIR}"/polap-bash-create-depth-file.sh 40 90 \
		"${_polap_var_mtcontigs_depth_range_graphfilter}"
	echo output: "${_polap_var_mtcontigs_depth_range_graphfilter}"
}

# Test: run-polap-r-depthfilter-gfa.R
# run-polap-function-seeds.sh
# function _polap_seeds_depthfilter-gfa
# function _run_polap_seeds-graph { # select seed contigs
test-r-depthfilter-gfa() {
	local _polap_var_mtcontigs=output

	# input
	local _polap_var_ga_contigger_edges_gfa=input/assembly_graph.gfa
	local _polap_var_mtcontigs_depth_range_graphfilter="${_polap_var_mtcontigs}/depth.range.graphfilter.txt"
	local _polap_var_mtcontigs_depth_range_graphfilter="input/1-depth.range.txt"
	echo input1: "${_polap_var_ga_contigger_edges_gfa}"
	echo input2: "${_polap_var_mtcontigs_depth_range_graphfilter}"

	# output
	local _polap_var_mtcontigs_gfa_all="${_polap_var_mtcontigs}/3-gfa.all.gfa"
	local _polap_var_mtcontigs_gfa_seq_filtered="${_polap_var_mtcontigs}/3-gfa.seq.depthfiltered.txt"
	local _polap_var_mtcontigs_gfa_seq_part="${_polap_var_mtcontigs}/3-gfa.seq.all.tsv"

	gfatools view \
		-S ${_polap_var_ga_contigger_edges_gfa} \
		>${_polap_var_mtcontigs_gfa_all}
	echo output1: ${_polap_var_mtcontigs_gfa_all}
	grep ^S ${_polap_var_mtcontigs_gfa_all} >${_polap_var_mtcontigs_gfa_seq_part}
	echo output2: ${_polap_var_mtcontigs_gfa_seq_part}
	Rscript ${_POLAPLIB_DIR}/run-polap-r-depthfilter-gfa.R \
		--gfa ${_polap_var_mtcontigs_gfa_seq_part} \
		--depth ${_polap_var_mtcontigs_depth_range_graphfilter} \
		--out ${_polap_var_mtcontigs_gfa_seq_filtered}
	echo output3: ${_polap_var_mtcontigs_gfa_seq_filtered}
	Rscript ${_POLAPLIB_DIR}/run-polap-r-depthfilter-gfa.R \
		--gfa ${_polap_var_mtcontigs_gfa_seq_part} \
		--lower-bound-depth 40 \
		--upper-bound-depth 90 \
		--out ${_polap_var_mtcontigs_gfa_seq_filtered}.copy
	echo Rscript ${_POLAPLIB_DIR}/run-polap-r-depthfilter-gfa.R \
		--gfa ${_polap_var_mtcontigs_gfa_seq_part} \
		--lower-bound-depth 40 \
		--upper-bound-depth 90 \
		--out ${_polap_var_mtcontigs_gfa_seq_filtered}.copy
	echo output4: ${_polap_var_mtcontigs_gfa_seq_filtered}.copy
}

test-polap-r-data-v2-alpha0() {
	# Plot for delta
	local _suppfigure_file="output/delta.pdf"
	Rscript ${_POLAPLIB_DIR}/run-polap-r-data-v2-alpha0.R \
		input/?.??.tsv \
		-l delta \
		-o "${_suppfigure_file}"
	echo open "${_suppfigure_file}"

	# Plot for alpha0
	local _suppfigure_file="output/alpha0.pdf"
	Rscript ${_POLAPLIB_DIR}/run-polap-r-data-v2-alpha0.R \
		input/?.??.tsv \
		-l alpha0 \
		-o "${_suppfigure_file}"
	echo open "${_suppfigure_file}"
}

test-polap-r-edges-stats() {
	local _polap_var_ga_contigger="input"
	local _polap_var_ga_contigger_edges_gfa="${_polap_var_ga_contigger}/graph_final.gfa"
	echo input1: "${_polap_var_ga_contigger_edges_gfa}"

	local _polap_var_ga_contigger="output"
	local _polap_var_ga_gfa_all="${_polap_var_ga_contigger}/3-gfa.all.gfa"
	local _polap_var_ga_gfa_seq_part="${_polap_var_ga_contigger}/3-gfa.seq.part.tsv"
	local _polap_var_ga_contigger_edges_stats="${_polap_var_ga_contigger}/edges_stats.txt"
	echo output1: "${_polap_var_ga_contigger_edges_stats}"

	gfatools view \
		-S ${_polap_var_ga_contigger_edges_gfa} \
		>${_polap_var_ga_gfa_all}
	grep ^S ${_polap_var_ga_gfa_all} >${_polap_var_ga_gfa_seq_part}
	Rscript ${_POLAPLIB_DIR}/run-polap-r-edges-stats.R \
		--gfa ${_polap_var_ga_gfa_seq_part} \
		--out ${_polap_var_ga_contigger_edges_stats}
}

test-polap-r-disassemble() {

	local A=(
		summary1-0.txt
		summary1-1.txt
		summary1-2.txt
		summary1-2-1.txt
		summary1-3.txt
	)

	for i in "${A[@]}"; do
		local j="${i%.txt}"
		Rscript ${_POLAPLIB_DIR}/run-polap-r-disassemble.R \
			--table input/$i \
			--out output/$j-ordered.txt \
			--plot input/$j-ordered.pdf
		echo output: "$j-ordered.txt"
	done
}

test-polap-r-preselect-annotation() {

	# a -> a437
	i=2
	Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-preselect-annotation.R \
		--table input/run-polap-r-preselect-annotation.input$i.txt --depth-range 27,27 --compare-mt-pt \
		--out output/run-polap-r-preselect-annotation.output$i.txt --gene-density 10 --plastid

	i=3
	Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-preselect-annotation.R \
		--table input/run-polap-r-preselect-annotation.input$i.txt --depth-range 37,37 --compare-mt-pt \
		--out output/run-polap-r-preselect-annotation.output$i.txt --gene-density 10 --plastid

	i=4
	Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-preselect-annotation.R \
		--table input/run-polap-r-preselect-annotation.input$i.txt --depth-range 51,312 --compare-mt-pt \
		--out output/run-polap-r-preselect-annotation.output$i.txt --gene-density 10 --plastid

	i=5
	Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-preselect-annotation.R \
		--table input/run-polap-r-preselect-annotation.input$i.txt --depth-range 58,384 --compare-mt-pt \
		--out output/run-polap-r-preselect-annotation.output$i.txt --gene-density 10 --plastid

	i=6
	Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-preselect-annotation.R \
		--table input/run-polap-r-preselect-annotation.input$i.txt --depth-range 73,420 --compare-mt-pt \
		--out output/run-polap-r-preselect-annotation.output$i.txt --gene-density 10 --plastid

	i=7
	Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-preselect-annotation.R \
		--table input/run-polap-r-preselect-annotation.input$i.txt --depth-range 78,486 --compare-mt-pt \
		--out output/run-polap-r-preselect-annotation.output$i.txt --gene-density 10 --plastid

	i=8
	Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-preselect-annotation.R \
		--table input/run-polap-r-preselect-annotation.input$i.txt --depth-range 111,555 --compare-mt-pt \
		--out output/run-polap-r-preselect-annotation.output$i.txt --gene-density 10 --plastid

	i=9
	Rscript --vanilla ${_POLAPLIB_DIR}/run-polap-r-preselect-annotation.R \
		--table input/run-polap-r-preselect-annotation.input$i.txt --depth-range 107,630 --compare-mt-pt \
		--out output/run-polap-r-preselect-annotation.output$i.txt --gene-density 10 --plastid

	# for i in {2..9}; do
	# 	cp a437/0/disassemble/1/1/$i/assembly_info_organelle_annotation_count-all.txt input/run-polap-r-preselect-annotation.input$i.txt
	# 	cp a437/0/disassemble/1/1/$i/51-mtcontigs/2/1-preselection.by.gene.density.txt output_expected/run-polap-r-preselect-annotation.output$i.txt
	# done
}

test-polap-py-find-cc-with-seeds() {
	for i in {2..9}; do
		# cp a/0/disassemble/1/1/$i/51-mtcontigs/2/4-gfa.depthfiltered.gfa input/run-polap-py-find-cc-with-seeds.input1.$i.gfa
		# cp a/0/disassemble/1/1/$i/51-mtcontigs/2/1-preselection.by.gene.density.txt input/run-polap-py-find-cc-with-seeds.input2.$i.txt
		# cp a/0/disassemble/1/1/$i/51-mtcontigs/2/5-gfa.depthfiltered.cc.seed.txt output_expected/run-polap-py-find-cc-with-seeds.output1.$i.txt
		python "${_POLAPLIB_DIR}"/polap-py-find-cc-with-seeds.py \
			--gfa input/run-polap-py-find-cc-with-seeds.input1.$i.gfa \
			--nodes input/run-polap-py-find-cc-with-seeds.input2.$i.txt \
			--output output/run-polap-py-find-cc-with-seeds.output1.$i.txt
	done
}

test-polap-py-find-plastid-gfa2fasta() {

	for i in {2..9}; do
		# cp -p a/0/disassemble/1/1/$i/30-contigger/graph_final.gfa input/run-polap-py-find-plastid-gfa2fasta.input1.$i.gfa
		# cp -p a/0/disassemble/1/1/$i/mt.contig.name input/run-polap-py-find-plastid-gfa2fasta.input2.$i.mt.contig.name
		# cp -pr a/0/disassemble/1/1/$i/52-mtdna output_expected/run-polap-py-find-plastid-gfa2fasta.output.$i.52-mtdna

		mkdir -p output/run-polap-py-find-plastid-gfa2fasta.output.$i.52-mtdna

		python "${_POLAPLIB_DIR}"/polap-py-find-plastid-gfa2fasta.py \
			--gfa input/run-polap-py-find-plastid-gfa2fasta.input1.$i.gfa \
			--seed input/run-polap-py-find-plastid-gfa2fasta.input2.$i.mt.contig.name \
			--out output/run-polap-py-find-plastid-gfa2fasta.output.$i.52-mtdna
	done

}

# example files needs to be changed because they point to something in a
# not input or output folder.
test-polap-py-unique-mtcontigs() {
	for i in {2..9}; do
		cp -p a/0/disassemble/1/1/$i/51-mtcontigs/file_hashes.tmp input/run-polap-py-unique-mtcontigs.input1.$i.tmp
		cp -p a/0/disassemble/1/1/$i/51-mtcontigs/filtered_files.tmp output_expected/run-polap-py-unique-mtcontigs.output1.$i.tmp
		mkdir -p output/run-polap-py-unique-mtcontigs.output1.$i
		python "${_POLAPLIB_DIR}"/polap-py-unique-mtcontigs.py \
			--input input/run-polap-py-unique-mtcontigs.input1.$i.tmp \
			--out output/run-polap-py-unique-mtcontigs.output1.$i.tmp \
			--prefix output/run-polap-py-unique-mtcontigs.output1.$i/mt.contig.name
	done
}

test-polap-py-compare2ptdna() {
	for i in {2..9}; do
		cp a2/0/disassemble/1/1/$i/52-mtdna/circular_path_1_concatenated.fa input/run-polap-py-compare2ptdna.input1.$i.fa
		cp -pr a2/0/disassemble/1/1/$i/52-mtdna/1 output_expected/run-polap-py-compare2ptdna.output1.$i.52-mtdna
		mkdir -p output/run-polap-py-compare2ptdna.output1.$i.52-mtdna/1
		python "${_POLAPLIB_DIR}"/polap-py-compare2ptdna.py \
			--seq1 input/ptdna-reference.fa \
			--seq2 input/run-polap-py-compare2ptdna.input1.$i.fa \
			--out output/run-polap-py-compare2ptdna.output1.$i.52-mtdna/1
		# a2/0/disassemble/1/1/$i/52-mtdna/1
	done
}

test-polap-r-mafft() {
	# bash "${_POLAPLIB_DIR}"/../polap.sh mafft-mtdna \
	# 	-a input/ptdna-reference.fa \
	# 	-b input/pt.subsample-polishing.reference.aligned.1.fa \
	# 	-o output/mafft1
	#
	# bash "${_POLAPLIB_DIR}"/../polap.sh mafft-mtdna \
	# 	-a input/ptdna-reference.fa \
	# 	-b input/pt.subsample-polishing.1.fa \
	# 	-o output/mafft2

	Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-mafft.R \
		--input output/mafft1/out.mafft \
		--out output/mafft1/out.txt

	Rscript --vanilla ${_POLAPLIB_DIR}/polap-r-mafft.R \
		--input output/mafft2/out.mafft \
		--out output/mafft2/out.txt
}

# test-bash-extract-three-edges-of-ptdna
# test-r-depth-distribution
# test-bash-create-depth-file
# test-r-depthfilter-gfa
# test-polap-r-data-v2-alpha0
# test-polap-r-edges-stats
# test-polap-r-disassemble
# test-polap-r-preselect-annotation
# test-polap-py-find-cc-with-seeds
# test-polap-py-find-plastid-gfa2fasta
# test-polap-py-unique-mtcontigs
# test-polap-py-compare2ptdna
test-polap-r-mafft
