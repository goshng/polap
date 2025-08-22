#!/usr/bin/awk -f
# Given sorted BED (B,Bp,C) and N (total length), output "As Ae" (1-based inclusive)
BEGIN{prev=0; best=0}
{
  if ($3 <= $2) next
  if ($2>prev) { s=prev+1; e=$2; L=e-s+1; if(L>best){best=L; As=s; Ae=e} }
  if ($3>prev) prev=$3
}
END{
  if (prev<N) { s=prev+1; e=N; L=e-s+1; if(L>best){best=L; As=s; Ae=e} }
  if (As==""||Ae=="") { As=1; Ae=1 }
  print As, Ae
}
