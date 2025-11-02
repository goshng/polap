# Input: CHROM  POS  DEPTH  (from `samtools depth -aa`)
# Output: contig  start  end  n_bases  sum_depth  mean_depth  sd_depth
BEGIN { OFS = "\t"; chr=""; bin=-1; B = (B?B:500); s=0; ss=0; n=0; start=1; }
NR==1 { chr=$1; bin=int(($2-1)/B); start=bin*B+1; }
{
  c=$1; p=$2; d=$3; b=int((p-1)/B);
  if (c != chr || b != bin) {
    end=(bin+1)*B;
    mean = n ? s/n : 0;
    var  = (n>1) ? (ss - s*s/n)/(n-1) : 0;
    sd   = (var>0)? sqrt(var): 0;
    print chr, start, end, n, s, mean, sd;
    # reset
    chr=c; bin=b; start=bin*B+1; s=0; ss=0; n=0;
  }
  s += d; ss += d*d; n++;
}
END {
  end=(bin+1)*B;
  mean = n ? s/n : 0;
  var  = (n>1) ? (ss - s*s/n)/(n-1) : 0;
  sd   = (var>0)? sqrt(var): 0;
  print chr, start, end, n, s, mean, sd;
}

