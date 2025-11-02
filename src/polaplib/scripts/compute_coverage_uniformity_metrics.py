#!/usr/bin/env python3
import sys, pandas as pd, numpy as np

if len(sys.argv) != 3:
    sys.stderr.write("Usage: compute_coverage_uniformity_metrics.py depth.tsv out.tsv\n")
    sys.exit(1)
inp, outp = sys.argv[1], sys.argv[2]

# depth.tsv: CHROM POS DEPTH
df = pd.read_csv(inp, sep='\t', header=None, names=['contig','pos','depth'], dtype={'contig':str,'pos':np.int64,'depth':np.float64})

def gini(x):
    x = np.asarray(x, dtype=np.float64)
    if x.size == 0: return np.nan
    if np.all(x == 0): return 0.0
    x_sorted = np.sort(x)
    n = x_sorted.size
    cumw = np.arange(1, n+1)
    G = (2.0 * np.sum(cumw * x_sorted)) / (n * x_sorted.sum()) - (n + 1.0) / n
    return float(G)

def summarize(group):
    x = group['depth'].to_numpy()
    n = x.size
    mu = float(np.mean(x)) if n else np.nan
    med = float(np.median(x)) if n else np.nan
    sd = float(np.std(x, ddof=1)) if n > 1 else 0.0
    mad = float(np.median(np.abs(x - med))) if n else np.nan
    cv = (sd / mu) if (mu and mu > 0) else np.nan
    q05, q20, q80, q95 = [float(np.quantile(x, q)) for q in (0.05,0.20,0.80,0.95)] if n else [np.nan]*4
    # Normalize by mean for fold metrics
    if mu and mu > 0:
        r = x / mu
        r20 = float(np.quantile(r, 0.20))
        fold80 = (1.0 / r20) if r20 > 0 else np.inf  # Picard definition
        pct_ge_0_5x = float(np.mean(x >= 0.5*mu))
        pct_ge_1x   = float(np.mean(x >= 1.0*mu))
        pct_ge_2x   = float(np.mean(x >= 2.0*mu))
        gg = gini(x)
    else:
        fold80 = np.nan; pct_ge_0_5x = pct_ge_1x = pct_ge_2x = np.nan; gg = np.nan
    return pd.Series(dict(
        n_bases=n, mean_depth=mu, median_depth=med, sd_depth=sd, mad_depth=mad, cv=cv,
        q05=q05, q20=q20, q80=q80, q95=q95, fold80=fold80,
        pct_ge_0_5x=pct_ge_0_5x, pct_ge_1x=pct_ge_1x, pct_ge_2x=pct_ge_2x, gini=gg
    ))

per = df.groupby('contig', sort=False).apply(summarize).reset_index()
# Global row
allrow = summarize(df)
allrow.name = 'ALL'
allrow = allrow.to_frame().T
allrow.insert(0, 'contig', 'ALL')

out = pd.concat([per, allrow], ignore_index=True)
out.to_csv(outp, sep='\t', index=False)

