#!/usr/bin/awk -f
# Version: v0.1.0
# prints id \t length for sequences in a FASTA (rough demo)
BEGIN{ id=""; len=0 }
$0 ~ /^>/ {
  if (id!="") print id "\t" len
  id=substr($0,2); sub(/[ \t].*$/,"",id)
  len=0; next
}
{ gsub(/[^A-Za-z]/,""); len+=length($0) }
END{ if (id!="") print id "\t" len }
