#!/usr/bin/env awk -f
#
# Filter PAF alignments by minimum aligned length and identity.
# Usage:
#   awk -v MINLEN=1200 -v MINID=0.88 -f filter_paf_by_identity_len.awk input.paf > reject_ids.txt
#
# PAF columns (0-based description, 1-based AWK fields):
#   $1  qname (query/read name)
#   $10 nmatch (# residue matches)
#   $11 aln_block_len (alignment block length)

BEGIN {
    if (MINLEN == "") MINLEN = 1200
    if (MINID  == "") MINID  = 0.88
    if (MINCOV == "") MINCOV = 0.30
}

# Skip comment lines or malformed lines
/^#/ || NF < 12 { next }

{
    nmatch = $10 + 0
    blen   = $11 + 0
    iden   = (blen > 0 ? nmatch / blen : 0)

    qlen=$2+0; qstart=$3+0; qend=$4+0
    cov=(qlen>0 ? (qend-qstart)/qlen : 0)

    if (blen >= MINLEN && iden >= MINID && cov >= MINCOV) {
        print $1
    }
}

#!/usr/bin/env awk -f
# Keep qname if: aln len >= MINLEN, identity >= MINID, and coverage >= MINCOV
# BEGIN{ if(MINLEN=="")MINLEN=1500; if(MINID=="")MINID=0.86; if(MINCOV=="")MINCOV=0.30 }
# NF<12{next}
# {
#   qlen=$2+0; qstart=$3+0; qend=$4+0
#   nmatch=$10+0; blen=$11+0; id=(blen>0?nmatch/blen:0)
#   cov=(qlen>0 ? (qend-qstart)/qlen : 0)
#
#   if(blen>=MINLEN && id>=MINID && cov>=MINCOV) print $1
# }

