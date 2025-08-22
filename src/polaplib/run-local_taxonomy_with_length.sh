#!/usr/bin/bash

# all_chloroplast_refseq.fa ptdna.fa

# run.sh chloroplast
# run.sh mitochondrion

type=$1

if [[ "$type" == "chloroplast" ]]; then
  infile=all_chloroplast_refseq.fa
  outfile=ptdna.fa
  if [[ ! -s "$infile" ]]; then
    esearch -db nucleotide -query 'chloroplast[Title] AND complete[Title] AND genome[Title] AND Viridiplantae[Organism] AND srcdb_refseq[PROP] NOT (clone[Title] OR cDNA[Title] OR mRNA[Title])' |
      efetch -format fasta >$infile
  fi
else
  infile=all_mitochondrion_refseq.fa
  outfile=mtdna.fa
  if [[ ! -s "$infile" ]]; then
    esearch -db nucleotide -query 'mitochondrion[Title] AND complete[Title] AND genome[Title] AND Viridiplantae[Organism] AND srcdb_refseq[PROP] NOT (clone[Title] OR cDNA[Title] OR mRNA[Title])' |
      efetch -format fasta >$infile
  fi
fi

echo "create a set of organelle sequences one per genus"
echo "input: $infile"
echo "output: $outfile"

bash local_taxonomy_with_length.sh $infile

Rscript local_taxonomy_with_length.R

cut -f1 taxonomy_deduplicated.tsv | tail -n +2 >selected_ids.txt
seqkit grep -f selected_ids.txt $infile >$outfile
