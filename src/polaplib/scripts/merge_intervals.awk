# merge_intervals.awk
# Input: sorted intervals "start<TAB>end"
NR==1 { S=$1; E=$2; next }
{
  if($1<=E){ if($2>E) E=$2 }
  else { print S"\t"E; S=$1; E=$2 }
}
END{ if(NR>0) print S"\t"E }
