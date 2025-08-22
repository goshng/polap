# percent_identity.awk
NR==1 {s1=$2; next}
NR==2 {s2=$2}
END {
  len=length(s1)
  if(len!=length(s2)){exit 1}
  for(i=1;i<=len;i++){
    a=substr(s1,i,1); b=substr(s2,i,1)
    A=toupper(a); B=toupper(b)
    if(a!="- " && b!="- "){ungap++; if(A==B)match++}
  }
  pid=(ungap>0)?100.0*match/ungap:0.0
  printf("%.6f\n", pid)
}