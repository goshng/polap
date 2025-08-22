#!/usr/bin/awk -f
# Build a single-read schematic if a read hits A,B,C,Bp with >=MIN overlap each
# Vars: sA,eA,sB,eB,sC,eC,sBp,eBp,MIN
BEGIN{OFS="\t"}
function ovl(a1,a2,b1,b2){ if (a1>b1) b1=a1; if (a2<b2) b2=a2; d=b2-b1; return (d>0?d:0) }
{
  q=$1; strand=$5; ts=$8+1; te=$9; if (ts>te){t=ts;ts=te;te=t}
  oa=ovl(ts,te,sA,eA); ob=ovl(ts,te,sB,eB); oc=ovl(ts,te,sC,eC); op=ovl(ts,te,sBp,eBp)
  if (oa>=MIN){ hit[q,"A"]=strand; sa[q]= (ts>sA?ts:sA); ea[q]= (te<eA?te:eA) }
  if (ob>=MIN){ hit[q,"B"]=strand; sb[q]= (ts>sB?ts:sB); eb[q]= (te<eB?te:eB) }
  if (oc>=MIN){ hit[q,"C"]=strand; sc[q]= (ts>sC?ts:sC); ec[q]= (te<eC?te:eC) }
  if (op>=MIN){ hit[q,"Bp"]=strand; sp[q]= (ts>sBp?ts:sBp); ep[q]= (te<eBp?te:eBp) }
  count[q]=0
}
END{
  chosen=""
  for (q in count){
    if ( (q SUBSEP "A") in hit && (q SUBSEP "B") in hit && (q SUBSEP "C") in hit && (q SUBSEP "Bp") in hit ){
      chosen=q; break
    }
  }
  if (chosen==""){ exit 1 }
  print "segment","start","end","strand","seg_start","seg_end"
  printf "A\t%d\t%d\t%s\t%d\t%d\n", sa[chosen], ea[chosen], hit[chosen,"A"], sA, eA
  printf "B\t%d\t%d\t%s\t%d\t%d\n", sb[chosen], eb[chosen], hit[chosen,"B"], sB, eB
  printf "C\t%d\t%d\t%s\t%d\t%d\n", sc[chosen], ec[chosen], hit[chosen,"C"], sC, eC
  printf "Bp\t%d\t%d\t%s\t%d\t%d\n", sp[chosen], ep[chosen], hit[chosen,"Bp"], sBp, eBp
}
