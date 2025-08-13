BEGIN { N=0 }
# only match header lines
/^>/ {
  print "edge_" ++N
}

