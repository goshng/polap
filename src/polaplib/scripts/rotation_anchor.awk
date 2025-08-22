# rotation_anchor.awk
# Input columns: qstart qend sstart send sstrand length pident qlen slen bitscore
function abs(x){return x<0?-x:x}
BEGIN{bestCoverLen=-1; bestNearDist=1e18; selS=1; haveCover=0}
{
  q1=$1+0; q2=$2+0; s1=$3+0; s2=$4+0; strand=$5; L=$6+0; qlen=$9+0; slen=$10+0
  qlo=(q1<q2?q1:q2); qhi=(q1>q2?q1:q2);
  if (qlo<=1 && 1<=qhi){
    off = 1 - q1;
    if (strand=="plus"){ s = s1 + off } else { s = s2 - off }
    if (s<1) s=1; if (s>slen) s=slen;
    if (L>bestCoverLen){bestCoverLen=L; selS=s; haveCover=1}
  } else if (!haveCover) {
    d = (abs(1-qlo) < abs(1-qhi)) ? abs(1-qlo) : abs(1-qhi);
    if (d < bestNearDist){
      if (abs(1-qlo) < abs(1-qhi)){
        s = (strand=="plus") ? s1 : s2
      } else {
        s = (strand=="plus") ? s2 : s1
      }
      if (s<1) s=1; if (s>slen) s=slen;
      bestNearDist=d; selS=s;
    }
  }
}
END{ print (selS>0?selS:1) }
