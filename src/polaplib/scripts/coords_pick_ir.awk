#!/usr/bin/awk -f
# coords_pick_ir.awk — pick longest reverse self-hit from show-coords -rclTH
BEGIN{FS="\t"; OFS="\t"; best=0}
NR==1 && $1 ~ /^[Ss]1$/ { next }       # skip header if present
$1 ~ /^[0-9]+$/ {
  s1=$1+0; e1=$2+0; s2=$3+0; e2=$4+0; len1=$5+0;
  # reverse orientation: E2 < S2 (and E1 > S1 typically). Require reverse.
  if (!(e2 < s2)) next;
  if (s1>e1){t=s1;s1=e1;e1=t}
  if (s2>e2){t=s2;s2=e2;e2=t}
  # skip full-diagonal self
  if (s1==s2 && e1==e2) next;
  L = (len1>0 ? len1 : (e1 - s1 + 1));
  if (L>best){best=L; bs1=s1; be1=e1; bs2=s2; be2=e2}
}
END{ if(best>0) printf "%d\t%d\t%d\t%d\n", bs1,be1,bs2,be2; else exit 1 }

