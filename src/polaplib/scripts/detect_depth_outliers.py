#!/usr/bin/env python3
import sys, pandas as pd, numpy as np

if len(sys.argv) != 3:
    sys.stderr.write("Usage: detect_depth_outliers.py bins.tsv outliers.bed\n"); sys.exit(1)
bins_tsv, out_bed = sys.argv[1], sys.argv[2]

# bins.tsv: contig start end n sum mean sd  (tab)
cols = ['contig','start','end','n','sum','mean','sd']
df = pd.read_csv(bins_tsv, sep='\t', header=None, names=cols, dtype={'contig':str})

out_rows = []
for contig, sub in df.groupby('contig', sort=False):
    m = float(np.median(sub['mean']))
    mad = float(np.median(np.abs(sub['mean'] - m)))
    if mad == 0:
        continue
    rz = 0.6745 * (sub['mean'] - m) / mad   # robust z
    sub = sub.assign(rz=rz, sign=np.where(rz>=0, 1, -1), is_out=np.abs(rz) >= 5.0)
    # merge adjacent outlier bins
    curr = None
    for row in sub.itertuples(index=False):
        if not row.is_out:
            if curr is not None:
                out_rows.append(curr); curr=None
            continue
        label = "HIGH" if row.sign>0 else "LOW"
        if curr is None:
            curr = dict(contig=row.contig, start=int(row.start), end=int(row.end), label=label,
                        n_bins=1, max_abs_z=abs(float(row.rz)))
        else:
            # extend if contiguous and same label
            if curr['end'] == int(row.start) and curr['label'] == label:
                curr['end'] = int(row.end)
                curr['n_bins'] += 1
                curr['max_abs_z'] = max(curr['max_abs_z'], abs(float(row.rz)))
            else:
                out_rows.append(curr)
                curr = dict(contig=row.contig, start=int(row.start), end=int(row.end),
                            label=label, n_bins=1, max_abs_z=abs(float(row.rz)))
    if curr is not None: out_rows.append(curr)

out = pd.DataFrame(out_rows, columns=['contig','start','end','label','n_bins','max_abs_z'])
out.to_csv(out_bed, sep='\t', header=False, index=False)

