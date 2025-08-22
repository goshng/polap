# compute_weighted_median.awk
# Usage: awk -v minlen=10000 -f compute_weighted_median.awk depth.summary.txt
# Input columns (mosdepth summary): chrom length bases mean min max
BEGIN { FS = OFS = "\t" }
NR==1 { next }  # skip header
$2 >= minlen { L[NR]=$2; D[NR]=$4; tot += $2 }
END {
  if (tot == 0) { print 0; exit }
  target = tot / 2.0
  cum = 0
  # file is assumed pre-sorted by mean depth ascending
  for (i=2; i<=NR; i++) {
    if (i in L) {
      cum += L[i]
      if (cum >= target) { print D[i]; exit }
    }
  }
  # fallback
  print D[NR]
}

