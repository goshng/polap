You're right — my classification picked LSC/SSC using `LN:i:`. If the actual `S`-line sequence lengths differ (or `LN` is missing/misleading), we must compute lengths from the **sequence itself**.

Here’s a drop‑in script that:

- reads `S` lines and prefers `len(sequence)` from the **sequence field**;
- if the sequence is `*`, it can optionally load real sequences from `--fasta` (plain FASTA; ids = segment names);
- only if neither is available, it falls back to `LN:i:`;
- expands `L` links with reverse‑complement twins;
- for cycles including the seed edge:

  - 1‑segment ring → `seed+`, `seed-`
  - 3‑segment → classifies IR (seen in both signs in any `L` line); among the other two, **LSC = longer by actual sequence length**, **SSC = shorter**; emits **LSC‑first** valid paths only.

```python
#!/usr/bin/env python3
import argparse, networkx as nx
from collections import defaultdict

REV = {"+":"-","-":"+"}
def base_id(n): return n[:-1]
def orient(n):  return n[-1]

def read_fasta(path):
    seqs = {}
    if not path: return seqs
    name=None; buf=[]
    with open(path) as fh:
        for ln in fh:
            if ln.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = ln[1:].strip().split()[0]
                buf = []
            else:
                buf.append(ln.strip())
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs

def parse_gfa(gfa_path, fasta_map=None, verbose=False):
    seg_len = defaultdict(int)      # final length preference: |sequence| > FASTA > LN:i:
    links   = set()
    G       = nx.DiGraph()
    with open(gfa_path) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"): continue
            p = ln.rstrip("\n").split("\t")
            if p[0] == "S":
                name = p[1]
                seq  = p[2]
                L_by_seq = len(seq) if seq != "*" else 0
                L_by_fa  = len(fasta_map.get(name, "")) if fasta_map else 0
                L_by_LN  = 0
                for tag in p[3:]:
                    if tag.startswith("LN:i:"):
                        try: L_by_LN = int(tag.split(":")[-1])
                        except: L_by_LN = 0
                        break
                # preference: sequence > fasta > LN
                L = L_by_seq or L_by_fa or L_by_LN
                seg_len[name] = L
                if verbose and L_by_seq==0 and L_by_fa==0 and L_by_LN>0:
                    print(f"[warn] {name}: using LN:i:{L_by_LN} (no sequence/FASTA)")
            elif p[0] == "L":
                f, fo, t, to = p[1], p[2], p[3], p[4]
                if fo in "+-" and to in "+-":
                    links.add((f,fo,t,to))
                    G.add_edge(f"{f}{fo}", f"{t}{to}")
                    links.add((t,REV[to],f,REV[fo]))
                    G.add_edge(f"{t}{REV[to]}", f"{f}{REV[fo]}")
    return seg_len, links, G

def oriented_cycles_including_seed(G, seed_base):
    return [cyc for cyc in nx.simple_cycles(G) if any(base_id(n)==seed_base for n in cyc)]

def compress_cycle_to_unique_bases(cyc):
    seen=set(); out=[]
    for n in cyc:
        b=base_id(n)
        if b not in seen:
            seen.add(b); out.append(n)
    return out

def orients_seen_in_links(links):
    seen = defaultdict(set)
    for f,fo,t,to in links:
        seen[f].add(fo); seen[t].add(to)
    return seen

def classify_three(bases, seg_len, orient_seen, verbose=False):
    # IR: appears in both '+' and '-' in any L rows; ties/none -> median by actual length
    cand = [b for b in bases if len(orient_seen[b]) >= 2]
    if len(cand) == 1:
        ir = cand[0]
    elif len(cand) > 1:
        lengths = {b: seg_len[b] for b in bases}
        med = sorted(lengths.values())[1]
        ir = min(cand, key=lambda b: abs(seg_len[b]-med))
    else:
        lengths = sorted([(seg_len[b], b) for b in bases])
        ir = lengths[1][1]
    others = [b for b in bases if b != ir]
    # Use actual sequence lengths to decide LSC/SSC
    if seg_len[others[0]] >= seg_len[others[1]]:
        lsc, ssc = others[0], others[1]
    else:
        lsc, ssc = others[1], others[0]
    if verbose:
        print(f"[classify] LSC={lsc}({seg_len[lsc]}), IR={ir}({seg_len[ir]}), SSC={ssc}({seg_len[ssc]})")
    return lsc, ir, ssc

def valid_paths_3seg(links, LSC, IR, SSC, verbose=False):
    out=[]
    for sA in "+-":
        for sB in "+-":
            if (LSC,sA,IR,sB) not in links:
                if verbose: print(f"[skip] missing L: {LSC}{sA}->{IR}{sB}")
                continue
            sBlast = REV[sB]
            ok=False
            for sC in "+-":
                if (SSC,sC,IR,sBlast) in links:
                    out.append([f"{LSC}{sA}", f"{IR}{sB}", f"{SSC}{sC}", f"{IR}{sBlast}"])
                    ok=True
                elif verbose:
                    print(f"[skip] missing L: {SSC}{sC}->{IR}{sBlast}")
            if not ok and verbose:
                print(f"[info] no SSC sign works for IR_last={IR}{sBlast}")
    # dedup
    seen=set(); uniq=[]
    for p in out:
        t=tuple(p)
        if t not in seen:
            uniq.append(p); seen.add(t)
    return uniq

def derive_paths(seg_len, links, G, seed_base, verbose=False):
    orient_seen = orients_seen_in_links(links)
    cycles = oriented_cycles_including_seed(G, seed_base)
    if verbose: print(f"[info] cycles including {seed_base}: {len(cycles)}")

    out=[]; seen=set()
    for i,cyc in enumerate(cycles,1):
        comp = compress_cycle_to_unique_bases(cyc)
        bases = []
        for n in comp:
            b=base_id(n)
            if b not in bases:
                bases.append(b)
        if verbose:
            print(f"[cycle#{i}] raw={cyc}")
            print(f"[cycle#{i}] uniq={comp} bases={bases}")

        if len(bases)==1 and bases[0]==seed_base:
            for s in "+-":
                tup=(f"{seed_base}{s}",)
                if tup not in seen:
                    out.append([f"{seed_base}{s}"]); seen.add(tup)
            continue

        if len(bases)!=3 or seed_base not in bases:
            if verbose: print(f"[cycle#{i}] ignored (nbases={len(bases)}, seed in? {seed_base in bases})")
            continue

        LSC, IR, SSC = classify_three(bases, seg_len, orient_seen, verbose=verbose)
        for p in valid_paths_3seg(links, LSC, IR, SSC, verbose=verbose):
            t=tuple(p)
            if t not in seen:
                out.append(p); seen.add(t)
    return out

def main():
    ap = argparse.ArgumentParser(description="Emit plastid circular paths including a given seed edge (1- or 3-seg). LSC/SSC decided by actual sequence length.")
    ap.add_argument("--gfa", required=True)
    ap.add_argument("--seed-edge", required=True, help="Edge ID (no sign)")
    ap.add_argument("--fasta", help="FASTA providing sequences for S lines that have '*'")
    ap.add_argument("--out", help="Output file (one path per line)")
    ap.add_argument("-v","--verbose", action="store_true")
    args = ap.parse_args()

    fa = read_fasta(args.fasta) if args.fasta else {}
    seg_len, links, G = parse_gfa(args.gfa, fasta_map=fa, verbose=args.verbose)
    paths = derive_paths(seg_len, links, G, args.seed_edge, verbose=args.verbose)

    if args.out:
        with open(args.out,"w") as fh:
            for p in paths:
                fh.write(",".join(p)+"\n")
    else:
        for p in paths:
            print(",".join(p))

    if args.verbose:
        print(f"[done] paths emitted: {len(paths)}")

if __name__ == "__main__":
    main()
```

### Why your output started with SSC

Your earlier run used `LN:i:` to order LSC vs SSC. If `edge_1`’s **actual sequence length** > `edge_2`’s `LN`, it could be mis‑classified as LSC, letting SSC‑first paths leak through. With the change above, we compare **real sequence lengths** from the `S` line (or from `--fasta` if `S` has `*`), so `edge_2` (LSC) will be the longer of the non‑IR pair and all emitted paths will start with `edge_2`.

### Usage

```bash
# When S lines contain sequences:
python plastid_seed_paths.py --gfa oatk.gfa --seed-edge edge_2 -v

# When S lines have '*' and you have a matching FASTA:
python plastid_seed_paths.py --gfa oatk.gfa --seed-edge edge_2 --fasta segments.fa -v
```

If you still see SSC‑first lines after this, run with `-v` and share the `[classify]` line — it will show the three lengths it used (LSC/IR/SSC) so we can spot any remaining mismatch.
