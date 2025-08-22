# flag_organelle.awk
# Usage: awk -v pivot=DEPTH -v mult=2.5 -v minlen=10000 -f flag_organelle.awk depth.summary.txt > summary.with_flag.tsv
BEGIN { FS = OFS = "\t" }
NR==1 { print $0, "org"; next }
{
  mean = $4; len = $2
  org = (len >= minlen && mean > mult * pivot) ? 1 : 0
  print $0, org
}

