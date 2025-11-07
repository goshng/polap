# FILE: scripts/paf_to_bed_with_pid.py
#!/usr/bin/env python3
# VERSION: 0.1.0
"""
From minimap2 PAF: keep hits by min_len, min_pid; emit BED with attributes on the target (mt) axis.

PID is computed preferentially as (1 - dv) if dv tag present, else (aln_len-NM)/aln_len if NM present,
else falls back to (1 - de), else NA (filtered out).
Refs: Li 2018 minimap2 (Bioinformatics).
"""
import sys, math


def parse_tags(fields):
    tags = {}
    for f in fields[12:]:
        try:
            k, t, v = f.split(":", 2)
            # dv/de float, NM int
            if t == "i":
                tags[k] = int(v)
            elif t in ("f", "Z"):
                try:
                    tags[k] = float(v)
                except:
                    tags[k] = v
            else:
                tags[k] = v
        except:
            pass
    return tags


if len(sys.argv) < 4:
    sys.stderr.write("Usage: paf_to_bed_with_pid.py <in.paf> <min_len> <min_pid>\n")
    sys.exit(1)

paf = sys.argv[1]
MIN_LEN = float(sys.argv[2])
MIN_PID = float(sys.argv[3])

with open(paf) as fh:
    for line in fh:
        if not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend = (
            f[0],
            int(f[1]),
            int(f[2]),
            int(f[3]),
            f[4],
            f[5],
            int(f[6]),
            int(f[7]),
            int(f[8]),
        )
        # alignment block length on target
        aln_len = abs(tend - tstart)
        tags = parse_tags(f)
        pid = None
        if "dv" in tags:
            pid = max(0.0, min(1.0, 1.0 - float(tags["dv"])))
        elif "NM" in tags and aln_len > 0:
            pid = max(0.0, min(1.0, (aln_len - float(tags["NM"])) / float(aln_len)))
        elif "de" in tags:
            pid = max(0.0, min(1.0, 1.0 - float(tags["de"])))
        else:
            continue
        if aln_len >= MIN_LEN and pid >= MIN_PID:
            # BED fields: chrom start end name score strand + attrs: pid aln_len qname qstart qend
            name = f"{qname}:{qstart}-{qend}"
            score = int(round(pid * 1000))
            print(
                "\t".join(
                    map(
                        str,
                        [
                            tname,
                            min(tstart, tend),
                            max(tstart, tend),
                            name,
                            score,
                            strand,
                            f"{pid:.4f}",
                            aln_len,
                            qname,
                            qstart,
                            qend,
                        ],
                    )
                )
            )
