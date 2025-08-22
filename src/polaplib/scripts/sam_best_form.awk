#!/usr/bin/awk -f
# From SAM (no header), pick best AS per read; count winners by rname(Form1..Form4)
# Var: MQ (min MAPQ)
BEGIN{FS="\t"; OFS="\t"}
{
  read=$1; rname=$3; mapq=$5+0; AS=0
  if (mapq < MQ) next
  for (i=12;i<=NF;i++) if ($i ~ /^AS:i:/){ split($i,a,":"); AS=a[3]+0 }
  if (AS==0) next
  if (!(read in best) || AS>best[read]) {best[read]=AS; which[read]=rname}
}
END{
  n1=n2=n3=n4=0
  for (r in which){
    if (which[r]=="Form1") n1++
    else if (which[r]=="Form2") n2++
    else if (which[r]=="Form3") n3++
    else if (which[r]=="Form4") n4++
  }
  tot12=n1+n2; tot34=n3+n4
  print "form","reads","fraction"
  printf "Form1\t%d\t%.6f\n", n1, (tot12? n1/tot12:0)
  printf "Form2\t%d\t%.6f\n", n2, (tot12? n2/tot12:0)
  printf "Form3\t%d\t%.6f\n", n3, (tot34? n3/tot34:0)
  printf "Form4\t%d\t%.6f\n", n4, (tot34? n4/tot34:0)
}
