# scripts/gfa_build_idmap.awk
# Version: v0.2.0
BEGIN{ FS=OFS="\t"; c=0 }
$1=="S" {
  old=$2
  if(!(old in seen)){ c++; seen[old]=c; print old, c }
}
END{
  if(c==0){ print "ERR: no S lines found" > "/dev/stderr"; exit 2 }
}
