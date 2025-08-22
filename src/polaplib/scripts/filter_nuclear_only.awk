# filter_nuclear_only.awk
# Usage: awk -v minlen=10000 -f filter_nuclear_only.awk summary.with_flag.tsv > nuclear.summary.tsv
BEGIN { FS = OFS = "\t" }
NR==1 { print; next }
$2 >= minlen && $NF == 0 { print }

