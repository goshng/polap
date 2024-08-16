#!/usr/bin/env bash

rm -rf o/50-annotation
mkdir -p o/50-annotation
grep -v '#' o/30-contigger/contigs_stats.txt | cut -f 1 >o/50-annotation/contig.name

# seqkit grep -f o/50-annotation/contig.name o/30-contigger/contigs.fasta -o o/50-annotation/contig.fasta
cp o/30-contigger/contigs.fasta o/50-annotation/contig.fasta
makeblastdb -dbtype nucl -in o/50-annotation/contig.fasta -out o/50-annotation/contig >/dev/null 2>&1

tblastn -query src/mt.1.c70.3.faa -db o/50-annotation/contig -out o/50-annotation/mtaa.blast -evalue 1e-30 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' -num_threads 56
src/run-polap-genes.R o/50-annotation/mtaa.blast o/50-annotation/mtaa.blast.bed >/dev/null 2>&1
sort -k1,1 -k2,2n o/50-annotation/mtaa.blast.bed >o/50-annotation/mtaa.blast.sorted.bed
mkdir o/50-annotation/mtaa.bed

while IFS= read -r contig; do
	grep -w $contig o/50-annotation/mtaa.blast.sorted.bed >o/50-annotation/mtaa.bed/"$contig".bed
	bedtools merge -i o/50-annotation/mtaa.bed/$contig.bed >o/50-annotation/mtaa.bed/"$contig".bed.txt
	printf "%s\t%d\n" $contig $(wc -l <o/50-annotation/mtaa.bed/$contig.bed.txt)
done <o/50-annotation/contig.name | sort -k2 -rn >o/50-annotation/mt.gene.count
tar zcf o/50-annotation/mtaa.bed.tar.gz o/50-annotation/mtaa.bed
rm -rf o/50-annotation/mtaa.bed

tblastn -query src/pt.2.c70.3.faa -db o/50-annotation/contig -out o/50-annotation/ptaa.blast -evalue 1e-30 -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles' -num_threads 56
src/run-polap-genes.R o/50-annotation/ptaa.blast o/50-annotation/ptaa.blast.bed >/dev/null 2>&1
sort -k1,1 -k2,2n o/50-annotation/ptaa.blast.bed >o/50-annotation/ptaa.blast.sorted.bed
mkdir o/50-annotation/ptaa.bed

while IFS= read -r contig; do
	grep -w $contig o/50-annotation/ptaa.blast.sorted.bed >o/50-annotation/ptaa.bed/"$contig".bed
	bedtools merge -i o/50-annotation/ptaa.bed/$contig.bed >o/50-annotation/ptaa.bed/"$contig".bed.txt
	printf "%s\t%d\n" $contig $(wc -l <o/50-annotation/ptaa.bed/$contig.bed.txt)
done <o/50-annotation/contig.name | sort -k2 -rn >o/50-annotation/pt.gene.count
tar zcf o/50-annotation/ptaa.bed.tar.gz o/50-annotation/ptaa.bed
rm -rf o/50-annotation/ptaa.bed

src/run-polap-mtcontig.R o \
	o/50-annotation/mt.contig.name \
	o/assembly_info_organelle_annotation_count.txt \
	--contigger \
	>/dev/null 2>&1
# o/assembly_info_organelle_annotation_count-all.txt \

echo "..."
echo "Step 3 has been finished!"
echo "  Select contigs that you suspect might have been from mitochondrial DNAs."
echo "  You can use:"
echo "    1. genome assembly graph"
echo "    2. presence of organelle gene"
echo "    3. copy numbers of contigs"
echo "  You may use the following for the three features of contigs."
echo "  Note that contigs are denoted by edge_<number>, e.g., edge_123"
echo "..."
echo "Whole-genome assembly graph: $PWD/o/30-contigger/graph_final.gfa"
echo "  See the assembly graph by using Bandage availabe at https://rrwick.github.io/Bandage/"
echo "Organelle annotation and copy numbers:"
echo "  $ column -t o/assembly_info_organelle_annotation_count.txt"
echo "edge_3" >o/50-annotation/mt.contig.name-1
echo "edge_1" >>o/50-annotation/mt.contig.name-1
echo "edge_2" >>o/50-annotation/mt.contig.name-1
echo "Step 4 is done manually by you."
echo "  See the table availabe at o/assembly_info_organelle_annotation_count.txt"
echo "  and the genome assembly graph."
echo "  You could watch YouTube: https://youtu.be/29OaiFCqAzI to select contigs."
echo "  The following is an example of the Flye assembly table with organelle gene counts."
echo "    V6 is the copy number for a contig."
echo "    mt is the number of mitochondrial genes."
echo "    pt is the number of plastid genes."
echo "    V9 is edeg numbers that constitute the contig at the corresponding line."
echo "$ column -t o/assembly_info_organelle_annotation_count.txt"
echo "V1        V2      V3  V4  V5  V6  V7    V8  V9         mt  pt"
echo "contig_4  239631  22  Y   N   1   none  *   4          33  2"
echo "contig_3  159195  22  N   N   1   none  *   3,1,-2,-1  56  83"
echo "  Edit o/50-annotation/mt.contig.name-1 to add lines of mtDNA contig candidates"
echo "edge_3"
echo "edge_1"
echo "edge_2"
echo "  You must have no empty lines in the mt.contig.name-1."
echo "Once you have your mt.contig.name file ready, then:"
echo Your next command is:
echo "$ src/run-polap-mtdna-carex-pseudochinensis-05.sh"
