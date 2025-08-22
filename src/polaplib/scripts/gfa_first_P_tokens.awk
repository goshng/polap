#!/usr/bin/awk -f
# Print the tokens field (col 3) from the first P-line
BEGIN{FS="\t"}
$1=="P"{ print $3; exit }
