# FILE: scripts/name_mtpts_from_hmmer.py
#!/usr/bin/env python3
# VERSION: 0.1.0
"""
Read tract TSV and hmmsearch --tblout; append ordered plastid gene list & tract name.

Input TSV columns (from merge_mt_hits_to_tracts.py):
tract_id  mt_contig  start  end  strand  n_hits  mean_pid  mean_alnlen  top_cp_qname  cp_min  cp_max

--tblout format: columns include target name (tract_id), query name (HMM/gene), hmm from-to (for order).

Output TSV: original columns + ordered_genes + name
Name = ordered genes joined by '-' if >=1 gene; else fallback cp[cp_min-cp_max].
"""
import sys, collections

if len(sys.argv) < 3:
    sys.stderr.write("Usage: name_mtpts_from_hmmer.py <tracts.tsv> <hmm.tblout>\n")
    sys.exit(1)

tsv = sys.argv[1]
tbl = sys.argv[2]

hits = collections.defaultdict(list)
with open(tbl) as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split()
        target = f[0]  # MTPT_xxx
        query = f[2]  # HMM/gene name
        try:
            qstart = int(f[17])
            qend = int(f[18])  # ali coord on target (hmmsearch --tblout spec)
        except:
            # fallback: collect without order
            qstart = 0
            qend = 0
        hits[target].append((min(qstart, qend), query))

print(
    "\t".join(
        [
            "tract_id",
            "mt_contig",
            "start",
            "end",
            "strand",
            "n_hits",
            "mean_pid",
            "mean_alnlen",
            "top_cp_qname",
            "cp_min",
            "cp_max",
            "ordered_genes",
            "name",
        ]
    )
)
with open(tsv) as fh:
    next(fh)  # header
    for line in fh:
        f = line.rstrip("\n").split("\t")
        tract_id = f[0]
        cpmin = f[9]
        cpmax = f[10]
        genelist = []
        if tract_id in hits and hits[tract_id]:
            ordered = sorted(hits[tract_id], key=lambda x: x[0])
            genelist = [g for _, g in ordered]
        ordered_str = ";".join(genelist) if genelist else ""
        if genelist:
            name = "-".join(genelist)
        else:
            name = f"cp[{f[8]}:{cpmin}-{cpmax}]"
        print("\t".join(f + [ordered_str, name]))
