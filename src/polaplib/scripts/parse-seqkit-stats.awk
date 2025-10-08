#!/usr/bin/awk -f
###############################################################################
# parse-seqkit-stats-v0.1.0.awk
# Version : v0.1.0
# Purpose : parse `seqkit stats -T -a` output safely
###############################################################################
BEGIN {
  FS = OFS = "\t"
}
NR == 1 {
  for (i=1; i<=NF; i++) {
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", $i)
    idx[$i] = i
  }
  next
}
NR >= 2 {
  tb  = (("sum_len" in idx)   ? $(idx["sum_len"])   : "NA")
  rc  = (("num_seqs" in idx)  ? $(idx["num_seqs"])  : "NA")
  ml  = (("avg_len" in idx)   ? $(idx["avg_len"])   : "NA")
  n50 = (("N50" in idx)       ? $(idx["N50"])       : "NA")
  mq  = (("AvgQual" in idx)   ? $(idx["AvgQual"])   : "NA")
  min = (("min_len" in idx)   ? $(idx["min_len"])   : "NA")
  max = (("max_len" in idx)   ? $(idx["max_len"])   : "NA")
  print sp, tb, rc, ml, n50, mq, min, max
  exit
}
