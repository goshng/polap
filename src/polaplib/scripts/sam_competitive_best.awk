#!/usr/bin/awk -f
# sam_competitive_best.awk
# Read SAM from stdin. For each read, keep the primary alignment with the highest AS.
# Output TSV: form  reads  fraction
BEGIN{
  FS = OFS = "\t";
  MQ = (MQ=="" ? 0 : MQ+0);  # default MQ=0 unless provided via -v MQ=...
}
$0 ~ /^@/ { next }  # skip headers
{
  q = $1; rname = $3; mapq = $5 + 0;
  if (rname == "*" || mapq < MQ) next;

  # extract AS:i: tag
  AS = 0;
  for (i=12; i<=NF; i++) {
    if ($i ~ /^AS:i:/) { split($i,a,":"); AS = a[3]+0; break; }
  }
  if (!(q in best) || AS > bestAS[q]) { best[q] = rname; bestAS[q] = AS; }
}
END{
  n1=n2=n3=n4=0;
  for (q in best) {
    f = best[q];
    if      (f=="Form1") n1++;
    else if (f=="Form2") n2++;
    else if (f=="Form3") n3++;
    else if (f=="Form4") n4++;
  }
  print "form","reads","fraction";
  tot12 = n1+n2; tot34 = n3+n4;
  print "Form1",n1,(tot12 ? n1/tot12 : 0);
  print "Form2",n2,(tot12 ? n2/tot12 : 0);
  print "Form3",n3,(tot34 ? n3/tot34 : 0);
  print "Form4",n4,(tot34 ? n4/tot34 : 0);
}

