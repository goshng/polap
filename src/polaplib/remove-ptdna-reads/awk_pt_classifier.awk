# awk_pt_classifier.awk  v0.1.0
# Input: PAF; Output: qname if ident>=MID, qcov>=QC, ALEN>=alen
BEGIN{FS=OFS="\t"}
{
  if (NF>=12) {
    al=$11+0; nm=$10+0; mq=$12+0
    id=(al>0? nm/al:0)
    qcv=($4-$3)/$2
    if (id>=MID && qcv>=QC && al>=ALEN) print $1
  }
}
