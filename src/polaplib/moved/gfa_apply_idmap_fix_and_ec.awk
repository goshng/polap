# scripts/gfa_apply_idmap_fix_and_ec.awk
# Version: v0.4.0
# Usage: awk -v EC_MODE=const|min|max|mean -v EC_CONST=1 -v EC_SCALE=1.0 -v EC_ROUND=round|ceil|floor \
#            -f gfa_apply_idmap_fix_and_ec.awk map.tsv input.gfa > normalized.gfa
#
# - Replaces IDs (S/L) using map.tsv.
# - On S: ensure ll:f using dp:i/dp:f if present else 1.0; also remember ll[new_id] for ec synthesis.
# - On L: overlap "*" => "0M"; ensure ec:i using selected EC_MODE from ll[u], ll[v].
#
BEGIN{
  FS=OFS="\t"
}
FNR==NR {
  idmap[$1]=$2
  next
}
# Helpers
function round_mode(x, mode,  r){
  if(mode=="ceil")   return int((x==int(x))?x:(x>0?int(x)+1:int(x)))
  else if(mode=="floor") return int((x==int(x))?x:(x<0?int(x)-1:int(x)))
  else { # round to nearest
    r = (x>=0)? int(x+0.5) : int(x-0.5)
    return r
  }
}
$1=="S" {
  old=$2
  if(!(old in idmap)){
    print "ERR: S id not in map:", old > "/dev/stderr"; exit 2
  }
  new=idmap[old]
  $2=new

  has_ll=0; has_dp=0; dpv="1.0"
  for(i=4;i<=NF;i++){
    if($i ~ /^ll:f:/) { has_ll=1; split($i,a,":"); ll[new]=a[3]+0 }
    else if($i ~ /^dp:i:/){ has_dp=1; split($i,a,":"); dpv=a[3] }
    else if($i ~ /^dp:f:/){ has_dp=1; split($i,a,":"); dpv=a[3] }
  }
  if(!has_ll){
    ll[new] = (has_dp ? dpv+0.0 : 1.0)
    $(++NF)="ll:f:" ll[new]
  }
  print; next
}
$1=="L" {
  # L <from> <fo> <to> <to> <ovl> [tags...]
  oldu=$2; oldv=$4
  if(!(oldu in idmap) || !(oldv in idmap)){
    print "ERR: L endpoint not in map:", oldu, oldv > "/dev/stderr"; exit 2
  }
  u=idmap[oldu]; v=idmap[oldv]
  $2=u; $4=v
  if($6=="*") $6="0M"

  # Build ec:i if missing
  has_ec=0
  for(i=7;i<=NF;i++){
    if($i ~ /^ec:i:/){ has_ec=1; break }
  }
  if(!has_ec){
    covu = (u in ll)? ll[u]:1.0
    covv = (v in ll)? ll[v]:1.0
    if(EC_MODE=="min")      ec = covu<covv?covu:covv
    else if(EC_MODE=="max") ec = covu>covv?covu:covv
    else if(EC_MODE=="mean") ec = (covu+covv)/2.0
    else { ec = EC_CONST+0.0 }  # const
    ec = ec * (EC_SCALE+0.0)
    ecint = round_mode(ec, EC_ROUND)
    if(ecint<0) ecint=0
    $(++NF) = "ec:i:" ecint
  }
  print; next
}
{ print }
