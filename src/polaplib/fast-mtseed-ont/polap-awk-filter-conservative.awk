# polap-awk-filter-conservative.awk  v0.0.1
# PAF: fields 0..11 (qname, qlen, qstart, qend, ..., nmatch, alen, mapq)
BEGIN{ FS=OFS="\t" }
{
  if (NF>=12) {
    alen = $11 + 0; nmatch = $10 + 0; mapq = $12 + 0
    ident = (alen>0? nmatch/alen : 0)
    qcov  = ($4 - $3) / ($2 + 0.0)
    if (ident>=0.85 && qcov>=0.20 && mapq>=10) print $0
  }
}
