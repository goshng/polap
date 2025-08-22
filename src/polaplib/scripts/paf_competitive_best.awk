#!/usr/bin/awk -f
# paf_competitive_best.awk
# Read PAF (minimap2 -c). Keep per-read best by aln length, then matches.
# Filters: MINLEN (alignment length) and MINID (matches/alnlen).
BEGIN{
  FS = OFS = "\t";
  if (MINLEN == "") MINLEN = 5000;
  if (MINID  == "") MINID  = 0.95;
}
{
  # PAF fields: qname(1) qlen(2) qstart(3) qend(4) strand(5) tname(6) tlen(7)
  #             tstart(8) tend(9) matches(10) alnlen(11) mapq(12) ...
  q = $1; t = $6;
  m = $10 + 0; al = $11 + 0;
  if (al < MINLEN) next;
  id = (al > 0 ? m/al : 0);
  if (id < MINID) next;

  if (!(q in best) || al > best_al[q] || (al == best_al[q] && m > best_m[q])) {
    best[q]    = t;
    best_al[q] = al;
    best_m[q]  = m;
  }
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

