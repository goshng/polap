# FILE: scripts/merge_mt_hits_to_tracts.py
#!/usr/bin/env python3
# VERSION: 0.1.0
"""
Merge overlapping cp->mt hits (BED with attrs) into tracts on the mt axis.
Outputs a TSV with rich attributes and prints to stdout.

TSV columns:
tract_id  mt_contig  start  end  strand  n_hits  mean_pid  mean_alnlen  top_cp_qname  cp_min  cp_max

Tract IDs are MTPT_<species>_<k>.
"""
import sys, collections

if len(sys.argv) < 4:
    sys.stderr.write(
        "Usage: merge_mt_hits_to_tracts.py <hits.bed> <merge_dist> <species>\n"
    )
    sys.exit(1)

bed = sys.argv[1]
merge_dist = int(sys.argv[2])
sp = sys.argv[3]


class Hit:
    __slots__ = (
        "t",
        "s",
        "e",
        "name",
        "score",
        "strand",
        "pid",
        "aln",
        "qname",
        "qs",
        "qe",
    )

    def __init__(self, row):
        self.t = row[0]
        self.s = int(row[1])
        self.e = int(row[2])
        self.name = row[3]
        self.score = int(row[4])
        self.strand = row[5]
        self.pid = float(row[6])
        self.aln = float(row[7])
        self.qname = row[8]
        self.qs = int(row[9])
        self.qe = int(row[10])


hits_by_contig = collections.defaultdict(list)
with open(bed) as fh:
    for line in fh:
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        hits_by_contig[f[0]].append(Hit(f))


def merge(hits):
    hits.sort(key=lambda h: (h.s, h.e))
    tr = []
    cur = [hits[0]]
    for h in hits[1:]:
        if h.s <= cur[-1].e + merge_dist:
            cur.append(h)
            if h.e > cur[-1].e:
                pass
        else:
            tr.append(cur)
            cur = [h]
    tr.append(cur)
    return tr


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
        ]
    )
)
k = 1
for t, arr in hits_by_contig.items():
    groups = merge(arr)
    for g in groups:
        start = min(h.s for h in g)
        end = max(h.e for h in g)
        # choose dominant strand (mode)
        strands = [h.strand for h in g]
        strand = max(set(strands), key=strands.count)
        mean_pid = sum(h.pid for h in g) / len(g)
        mean_aln = sum(h.aln for h in g) / len(g)
        cp_names = [h.qname for h in g]
        top_cp = max(set(cp_names), key=cp_names.count)
        cp_min = min(min(h.qs, h.qe) for h in g)
        cp_max = max(max(h.qs, h.qe) for h in g)
        tract_id = f"MTPT_{sp}_{k:03d}"
        print(
            "\t".join(
                map(
                    str,
                    [
                        tract_id,
                        t,
                        start,
                        end,
                        strand,
                        len(g),
                        f"{mean_pid:.4f}",
                        f"{mean_aln:.1f}",
                        top_cp,
                        cp_min,
                        cp_max,
                    ],
                )
            )
        )
        k += 1
