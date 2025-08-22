#!/usr/bin/awk -f
# Extract S-lines from GFA as FASTA with the segment name as header
BEGIN{FS="\t"}
$1=="S"{
  print ">" $2
  print $3
}
