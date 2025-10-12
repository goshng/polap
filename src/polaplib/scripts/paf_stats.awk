# Version: v0.5.1
# Usage: awk -f paf_stats.awk in.paf > out.tsv
BEGIN{ FS=OFS="\t"; print "n_align","mean_identity","total_aligned" }
{
  nmatch=$10+0; alen=$11+0;
  if (alen>0){ sum_id+=(nmatch/alen) }
  sum_al+=alen; n++
}
END{
  mi = (n>0)? sum_id/n : 0.0;
  printf "%d\t%.6f\t%d\n", n, mi, sum_al
}
