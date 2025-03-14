flye_root=/home/goshng/all/polap/Flye
# flye/tests/data/ecoli_500kb_reads.fastq.gz
# reads_file=flye/tests/data/ecoli_500kb.fasta
# reads_file=flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz

out_dir=${flye_root}/mt-d
mt_fastq=/home/goshng/all/polap/aflye1/Macadamia_tetraphylla/o/2/03-seeds/ptgaul/0.fq.gz
mt_fastq=${flye_root}/mt/1.fq.gz
mt_fastq=/home/goshng/all/polap/aflye1/Macadamia_tetraphylla/o/2/04-subsample/ptgaul/3.fq.gz
mt_fastq=/home/goshng/all/polap/aflye1/Punica_granatum/o/3/04-subsample/ptgaul/0.fq.gz
mt_fastq=/home/goshng/all/polap/aflye1/Lolium_perenne/o/2/04-subsample/ptgaul/0.fq.gz
mt_fastq=/home/goshng/all/polap/aflye1/Lolium_perenne/o/2/03-seeds/ptgaul/12.fq.gz
rm -rf ${out_dir}
source $HOME/miniconda3/bin/activate polap
bin/dflye \
	--nano-raw ${mt_fastq} \
	--debug \
	--stop-after contigger \
	--asm-coverage 30 \
	-g 800000 -o $out_dir -t 56 -m 10000
exit
bin/dflye \
	--directional-reads

--stop-after contigger

flye_root=/home/goshng/all/polap/Flye
out_dir=${flye_root}/o-d
ecoli_hifi_fastq=${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz
hifi_fastq=${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz
both_fastq=${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz
forward_fastq=${out_dir}/00-assembly/forward.fq
forward_paf=${out_dir}/00-assembly/forward.paf

# repeat graph with the forward filtered data
out_dir=${flye_root}/mt
disjointigs=${out_dir}/c0.fa
both_fastq=${out_dir}/lk.fq.gz
rm -rf ${out_dir}/20-repeat
mkdir -p ${out_dir}/20-repeat
bin/dflye-modules repeat \
	--disjointigs ${disjointigs} \
	--reads ${both_fastq} \
	--out-dir ${out_dir}/20-repeat \
	--config ${flye_root}/flye/config/bin_cfg/asm_raw_reads.cfg \
	--log ${out_dir}/flye.log \
	--threads 8 \
	--debug \
	--min-ovlp 3000 \
	--directional-reads
exit

rm -rf ${out_dir}
bin/dflye \
	--pacbio-corr ${hifi_fastq} \
	--debug \
	--directional-reads \
	-g 500k -o $out_dir -t 8 -m 1000
exit

# repeat graph with the forward filtered data
rm -rf ${out_dir}/20-repeat
mkdir -p ${out_dir}/20-repeat
bin/dflye-modules repeat \
	--disjointigs ${out_dir}/path1.fasta \
	--reads ${both_fastq} \
	--out-dir ${out_dir}/20-repeat \
	--config ${flye_root}/flye/config/bin_cfg/asm_corrected_reads.cfg \
	--log ${out_dir}/flye.log \
	--threads 8 \
	--debug \
	--min-ovlp 1000 \
	--directional-reads
exit

# Filter out reverse-complementary reads
bin/dflye-minimap2 \
	${out_dir}/path1.fasta \
	${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz \
	-x map-pb -t 8 -a -p 0.5 -N 10 \
	--sam-hit-only -L -K 1G -z 1000 -Q \
	--secondary-seq -I 64G |
	bin/dflye-samtools fastq -G 16 - >${forward_fastq}
exit

# Check the forward mapping
bin/dflye-minimap2 \
	${out_dir}/path1.fasta \
	${forward_fastq} \
	-x map-pb -t 8 -p 0.5 -N 10 \
	--sam-hit-only -L -K 1G -z 1000 -Q \
	--secondary-seq -I 64G >${forward_paf}
exit

# All
#
rm -rf ${out_dir}
bin/dflye \
	--pacbio-corr ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz \
	--stop-after consensus \
	--directional-reads \
	--debug \
	-g 500k -o $out_dir -t 8 -m 1000
exit

flye_root=/home/goshng/all/polap/Flye
out_dir=${flye_root}/o-d
rm -rf ${out_dir}
mkdir -p \
	${out_dir}/00-assembly \
	${out_dir}/10-consensus \
	${out_dir}/20-repeat \
	${out_dir}/30-contigger \
	${out_dir}/40-polishing

# other data
flye_root=/home/goshng/all/polap/Flye
out_dir=${flye_root}/l-d
mkdir -p \
	${out_dir}/00-assembly \
	${out_dir}/10-consensus \
	${out_dir}/20-repeat \
	${out_dir}/30-contigger \
	${out_dir}/40-polishing

bin/dflye-modules assemble \
	--reads ${flye_root}/l.fq \
	--out-asm ${out_dir}/00-assembly/draft_assembly.fasta \
	--config ${flye_root}/flye/config/bin_cfg/asm_raw_reads.cfg \
	--log ${out_dir}/flye.log \
	--threads 8 \
	--debug \
	--genome-size 500000 \
	--min-ovlp 1000 \
	--directional-reads

exit

# Stage: assembly
# input1: ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz
# output: ${out_dir}/00-assembly/draft_assembly.fasta
#
bin/dflye-modules assemble \
	--reads ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz \
	--out-asm ${out_dir}/00-assembly/draft_assembly.fasta \
	--config ${flye_root}/flye/config/bin_cfg/asm_corrected_reads.cfg \
	--log ${out_dir}/flye.log \
	--threads 8 \
	--debug \
	--genome-size 500000 \
	--min-ovlp 1000 \
	--directional-reads

# Stage: consensus 1
# input: ${out_dir}/00-assembly/draft_assembly.fasta
# output: ${out_dir}/10-consensus/minimap.bam
#
bin/dflye-minimap2 \
	${out_dir}/00-assembly/draft_assembly.fasta \
	${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz \
	-x map-pb -t 8 -a -p 0.5 -N 10 \
	--sam-hit-only -L -K 1G -z 1000 -Q \
	--secondary-seq -I 64G |
	bin/dflye-samtools view -T ${out_dir}/00-assembly/draft_assembly.fasta -u - |
	bin/dflye-samtools view -h - |                             # dflye:
	awk '{if ($1 ~ /^@/ || int($2 / 16) % 2 == 0) print $0}' | # dflye: Filter for forward-forward reads
	bin/dflye-samtools view -bS - |                            # dflye: Convert filtered SAM back to BAM
	bin/dflye-samtools sort -T ${out_dir}/10-consensus/sort_x_y \
		-O bam -@ 4 -l 1 -m 1G -o ${out_dir}/10-consensus/minimap.bam
exit

# Stage: consensus 2
# Compute consensus -> hard
# minimap.bam -> consensus.fasta
# input1: ${out_dir}/10-consensus/minimap.bam
# input1: ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz
# output: ${out_dir}/10-consensus/consensus.fasta

# Filter out reverse-complementary reads
bin/dflye-minimap2 \
	${out_dir}/00-assembly/draft_assembly.fasta \
	${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz \
	-x map-pb -t 8 -a -p 0.5 -N 10 \
	--sam-hit-only -L -K 1G -z 1000 -Q \
	--secondary-seq -I 64G |
	bin/dflye-samtools fastq -G 16 - >${out_dir}/00-assembly/forward.fq
exit

# Stage: repeat
# input1: ${out_dir}/10-consensus/consensus.fasta
# output: ${out_dir}/20-repeat
#
bin/dflye-modules repeat \
	--disjointigs ${out_dir}/10-consensus/consensus.fasta \
	--reads ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz \
	--out-dir ${out_dir}/20-repeat \
	--config ${flye_root}/flye/config/bin_cfg/asm_corrected_reads.cfg \
	--log ${out_dir}/flye.log \
	--threads 8 \
	--debug \
	--min-ovlp 1000 \
	--directional-reads

# Stage: contigger
# input1: ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz
# input2: ${flye_root}/flye/config/bin_cfg/asm_corrected_reads.cfg
# input3: ${out_dir}/20-repeat/repeat_graph_edges.fasta
# input4: ${out_dir}/20-repeat/repeat_graph_dump
# input5: ${out_dir}/20-repeat/read_alignment_dump
# output: ${out_dir}/30-contigger
#
mkdir -p ${out_dir}/30-contigger
bin/dflye-modules contigger \
	--graph-edges ${out_dir}/20-repeat/repeat_graph_edges.fasta \
	--reads ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz \
	--out-dir ${out_dir}/30-contigger \
	--config ${flye_root}/flye/config/bin_cfg/asm_corrected_reads.cfg \
	--repeat-graph ${out_dir}/20-repeat/repeat_graph_dump \
	--graph-aln ${out_dir}/20-repeat/read_alignment_dump \
	--log ${out_dir}/flye.log \
	--threads 8 \
	--debug \
	--min-ovlp 1000 \
	--directional-reads
exit

# Stage: polishing 1
# input1: ${out_dir}/30-contigger/contigs.fasta
# input2: ${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz
# output: ${out_dir}/40-polishing/minimap_1.bam
#
bin/dflye-minimap2 \
	${out_dir}/30-contigger/contigs.fasta \
	${flye_root}/flye/tests/data/ecoli_500kb_reads_hifi.fastq.gz -x map-pb -t 8 -a -p 0.5 -N 10 \
	--sam-hit-only -L -K 1G -z 1000 -Q \
	--secondary-seq -I 64G |
	flye-samtools view -T ${out_dir}/30-contigger/contigs.fasta -u - |
	flye-samtools sort -T ${out_dir}/40-polishing/sort_241207_112408 -O bam -@ 4 -l 1 -m 1G \
		-o ${out_dir}/40-polishing/minimap_1.bam

# Stage: polishing 2
# input1: ${out_dir}/40-polishing/filtered_contigs.fasta
# input1: ${out_dir}/30-contigger/graph_final.fasta
# output: ${out_dir}/40-polishing/edges_aln.bam
#
bin/dflye-minimap2 \
	${out_dir}/40-polishing/filtered_contigs.fasta \
	${out_dir}/30-contigger/graph_final.fasta -x map-pb -t 8 -a -p 0.5 -N 10 \
	--sam-hit-only -L -K 1G -z 1000 -Q \
	--secondary-seq -I 64G |
	flye-samtools view -T ${out_dir}/40-polishing/filtered_contigs.fasta -u - |
	flye-samtools sort -T ${out_dir}/40-polishing/sort_241207_112414 -O bam -@ 4 -l 1 -m 1G \
		-o ${out_dir}/40-polishing/edges_aln.bam

# Stage: finalize
#
