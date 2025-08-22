# extract_seq.awk
/^S/ { print ">"$2"\n"$3 }