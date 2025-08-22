#!/usr/bin/awk -f
# Select the longest reverse-strand self-hit with non-overlapping intervals
# Inputs: PAF; Vars: ML (min length, default 8000), MI (min identity %, default 0)
BEGIN{OFS="\t"; if(ML=="") ML=8000; if(MI=="") MI=0; best=0}
$5=="-" && $1==$6 {
  alen = ($11+0>0)?$11:$10
  nm = $10+0
  id = (alen>0)? (100.0*nm/alen) : 0
  if (alen>=ML && id>=MI) {
    s1=$8+1; e1=$9; if(s1>e1){t=s1;s1=e1;e1=t}
    s2=$3+1; e2=$4; if(s2>e2){t=s2;s2=e2;e2=t}
    if (e1<s2 || e2<s1) {
      if (alen>best){best=alen; bs1=s1; be1=e1; bs2=s2; be2=e2}
    }
  }
}
END{ if(best>0) printf "%d\t%d\t%d\t%d\n", bs1,be1,bs2,be2 }
