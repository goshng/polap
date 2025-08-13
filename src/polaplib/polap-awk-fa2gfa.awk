BEGIN {
  print "H\tVN:Z:1.0"
}
function flush_record(   edge, path) {
  if (!in_rec) return
  edge = "edge_" ++edge_i
  S_lines[edge_i] = "S\t" edge "\t" seq "\tdp:i:1"
  path = "contig_" edge_i
  P_lines[edge_i] = "P\t" path "\t" edge "+\t*"
  seq = ""
  in_rec = 0
}
# header line
/^>/ {
  flush_record()
  in_rec = 1
  next
}
# sequence lines (multi-line OK)
{
  gsub(/\r/, "")
  seq = seq $0
}
END {
  flush_record()
  for (i = 1; i <= edge_i; i++) print S_lines[i]
  for (i = 1; i <= edge_i; i++) print P_lines[i]
}
