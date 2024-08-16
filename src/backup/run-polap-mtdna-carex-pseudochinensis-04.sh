#!/usr/bin/env bash

rm -rf o/50-annotation
mkdir -p o/50-annotation

if [ ! -s mitochondrion.1.protein.faa ]; then
	wget https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.1.protein.faa.gz
	wget https://ftp.ncbi.nlm.nih.gov/refseq/release/plastid/plastid.2.protein.faa.gz
	gunzip *.gz
fi

cp $PWD/o/30-contigger/contigs.fasta o/50-annotation
makeblastdb -in o/50-annotation/contigs.fasta -dbtype nucl -out contigs
tblastn -query mitochondrion.1.protein.faa -db contigs -out mtaa.blast -evalue 1e-30 -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles" -num_threads 56
tblastn -query plastid.2.protein.faa -db contigs -out ptaa.blast -evalue 1e-30 -outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles" -num_threads 56

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
echo "  $ column -t $PWD/o/30-contigger/contigs_stats.txt"
echo "edge_3" >o/50-annotation/mt.contig.name-1
echo "edge_1" >>o/50-annotation/mt.contig.name-1
echo "edge_2" >>o/50-annotation/mt.contig.name-1
echo "Step 4 is done manually by you."
echo "..."
echo "library(readr)"
echo "library(dplyr)"
echo "x <- read_tsv(\"mtaa.blast\", col_names=F)"
echo "x %>% count(X2) %>% arrange(desc(n)) %>% print(n=30)"
echo "..."
echo "  See the table availabe at $PWD/o/30-contigger/contigs_stats.txt"
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
echo "$ src/run-polap-mtdna-carex-pseudochinensis-05.sh <MT.CONTIG.NAME-NUMBER>"
echo MT.CONTIG.NAME-NUMBER is an integer: e.g., 1.
