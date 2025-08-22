# plus_stats.awk
# MODE=all  -> prints one number: sum_length_all
# MODE=plus -> prints: "<sum_length> <weighted_identity> <best_bitscore>"
# Input columns: qstart qend sstart send sstrand length pident qlen slen bitscore
BEGIN{ s=0; w=0; best=0; mode_all=(MODE=="all")?1:0 }
{
  if (mode_all) { s += $6 }
  else { s += $6; w += $6*$7; if ($10>best) best=$10 }
}
END{
  if (mode_all) { printf("%d\n", s+0) }
  else {
    pid = (s>0)?(w/s):0.0
    printf("%.0f %.6f %.0f\n", s, pid, best)
  }
}
