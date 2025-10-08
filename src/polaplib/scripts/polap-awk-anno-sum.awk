#!/usr/bin/awk -f
# polap-awk-anno-sum.awk
# Version : v1.1.0  (2025-10-07)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Input  : contig-annotation-depth-table (space/tab separated)
# Header must include: Length, MT, PT (order agnostic).
# Output : a single TSV line:
#          total_len_bp  total_len_kb  mt_genes_sum  pt_genes_sum  n_contigs

BEGIN { FS="[ \t]+"; OFS="\t"; L=MT=PT=0; lsum=mt=pt=n=0 }
NR==1{
  for(i=1;i<=NF;i++){
    if($i=="Length") L=i;
    else if($i=="MT") MT=i;
    else if($i=="PT") PT=i;
  }
  if(!L){ print "[ERR] Missing Length column" > "/dev/stderr"; exit 2 }
  next
}
NR>1{
  if($0 ~ /^[[:space:]]*$/) next
  lsum += (L ? $(L)+0 : 0)
  mt   += (MT? $(MT)+0: 0)
  pt   += (PT? $(PT)+0: 0)
  n++
}
END{
  printf "%.0f\t%.0f\t%d\t%d\t%d\n", lsum, (lsum>0?lsum/1000.0:0), mt, pt, n
}
