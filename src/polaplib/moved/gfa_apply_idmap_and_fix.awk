# scripts/gfa_apply_idmap_and_fix.awk
# Version: v0.3.0
# Usage: awk -f gfa_apply_idmap_and_fix.awk map.tsv input.gfa > normalized.gfa
# Behavior:
#   - Replace S/L endpoint IDs using map.tsv (old -> new numeric).
#   - On L, if overlap == "*", set to "0M".
#   - On S, add ll:f from dp:i or dp:f if ll:f missing; default ll:f:1.0 otherwise.

BEGIN{
  FS=OFS="\t"
  # First file is the map
  reading_map=1
}
FNR==NR {
  idmap[$1]=$2
  next
}
{
  # Now processing the GFA
  if($1=="S"){
    name=$2
    if(!(name in idmap)){
      print "ERR: S id not in map:", name > "/dev/stderr"; exit 2
    }
    $2=idmap[name]
    # ll:f tag handling
    has_ll=0; has_dp=0; dpv="1.0"
    for(i=4;i<=NF;i++){
      if($i ~ /^ll:f:/) has_ll=1
      else if($i ~ /^dp:i:/){ split($i,a,":"); has_dp=1; dpv=a[3] }
      else if($i ~ /^dp:f:/){ split($i,a,":"); has_dp=1; dpv=a[3] }
    }
    if(!has_ll){
      if(has_dp) $(++NF)="ll:f:" dpv
      else       $(++NF)="ll:f:1.0"
    }
    print; next
  }
  else if($1=="L"){
    # L <from> <f_or> <to> <t_or> <ovl> ...
    u=$2; v=$4
    if(!(u in idmap) || !(v in idmap)){
      print "ERR: L endpoint not in map:", u, v > "/dev/stderr"; exit 2
    }
    $2=idmap[u]
    $4=idmap[v]
    if($6=="*") $6="0M"
    print; next
  }
  else{
    print; next
  }
}
