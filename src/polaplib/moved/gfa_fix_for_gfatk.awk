# scripts/gfa_fix_for_gfatk.awk
# Version: v0.1.0
# Goals:
#  - Ensure GFA1 H line exists and pass-through.
#  - On S lines: add ll:f from dp:i/f if ll:f missing.
#  - On L lines: if Overlap field == "*", replace with "0M".
#    (gfatk expects simple <int>M overlaps.)
BEGIN { FS=OFS="\t" }
$1=="H" { print; next }
$1=="S" {
  # Fields: S <name> <seq> [tags...]
  # Add ll:f if missing; derive from dp:i or dp:f when present.
  has_ll = 0; has_dp_i=0; has_dp_f=0; dpv=""
  for(i=4; i<=NF; i++){
    if ($i ~ /^ll:f:/) has_ll=1
    else if ($i ~ /^dp:i:/) { has_dp_i=1; split($i,a,":"); dpv=a[3] }
    else if ($i ~ /^dp:f:/) { has_dp_f=1; split($i,a,":"); dpv=a[3] }
  }
  if (!has_ll) {
    if (has_dp_i || has_dp_f) $(++NF)="ll:f:" dpv
    else                      $(++NF)="ll:f:1.0"
  }
  print; next
}
$1=="L" {
  # Fields: L <from> <orient> <to> <orient> <overlap> [tags...]
  if ($6=="*") $6="0M"
  print; next
}
{ print }
