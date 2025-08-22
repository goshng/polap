# sum_len_times_mean.awk
# Usage: awk -f sum_len_times_mean.awk depth.summary.txt
BEGIN { FS = OFS = "\t" }
NR>1 { sum += ($2 * $4) }
END { print sum+0 }

