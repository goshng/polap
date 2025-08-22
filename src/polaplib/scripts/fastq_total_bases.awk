# fastq_total_bases.awk
# Usage: awk -f fastq_total_bases.awk reads.fastq > total_bases.txt
(NR%4)==2 { L += length($0) }
END { print L+0 }

