#!/usr/bin/env python3
import sys, argparse, gzip, math


def open_text(fn):
    if fn == "-" or fn is None:
        return sys.stdin
    return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)


def rc(seq):
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def read_fx(fn):
    # Simple FASTA/FASTQ reader, returns dict[name]=sequence
    seqs = {}
    with open_text(fn) as f:
        line = f.readline()
        if not line:
            return seqs
        is_fastq = line.startswith("@")
        hdr = None
        buf = []

        def commit():
            nonlocal hdr, buf
            if hdr is None:
                return
            name = hdr.split()[0][1:]
            if name not in seqs:
                seqs[name] = "".join(buf)
            hdr = None
            buf = []

        if line[0] in ">@":
            hdr = line
        else:
            raise ValueError("Input must be FASTA or FASTQ")
        for line in f:
            c = line[:1]
            if c in ">@":
                commit()
                hdr = line
                buf = []
            elif is_fastq and c == "+":
                # skip qualities
                seqlen = sum(len(x.strip()) for x in buf)
                got = 0
                for qline in f:
                    got += len(qline.strip())
                    if got >= seqlen:
                        break
                commit()
                hdr = None
                # read next header if any
                nxt = f.readline()
                if not nxt:
                    break
                if nxt[:1] in "@>":
                    hdr = nxt
                    buf = []
                else:
                    raise ValueError("FASTQ parsing error")
            else:
                if c != "\n":
                    buf.append(line.strip())
        commit()
    return seqs


def approx_identity(nmatch, alen):
    try:
        return float(nmatch) / float(alen)
    except ZeroDivisionError:
        return 0.0


def parse_paf(paf_path, min_id, min_ovl):
    # Return edges: for each read, a list of candidate extensions (to, strand, ovl, iden)
    edges = {}
    with open_text(paf_path) as f:
        for ln in f:
            if not ln or ln[0] == "#":
                continue
            a = ln.rstrip("\n").split("\t")
            if len(a) < 12:
                continue
            qn, ql, qs, qe, st, tn, tl, ts, te, nm, al, mapq = a[:12]
            ql = int(ql)
            qs = int(qs)
            qe = int(qe)
            tl = int(tl)
            ts = int(ts)
            te = int(te)
            nm = int(nm)
            al = int(al)
            if qn == tn:
                continue
            ovl = min(qe - qs, te - ts)
            if ovl < min_ovl:
                continue
            iden = approx_identity(nm, al)
            if iden < min_id:
                continue
            edges.setdefault(qn, []).append((tn, st, ovl, iden))
            # also useful to know reverse direction when we extend from tn back into qn
            edges.setdefault(tn, []).append((qn, "+" if st == "+" else "-", ovl, iden))
    # sort by quality (iden * ovl)
    for k in edges:
        edges[k].sort(key=lambda x: (x[3], x[2]), reverse=True)
    return edges


def choose_next(current, used, edges):
    # pick the best neighbor not used yet
    for tn, st, ovl, iden in edges.get(current, []):
        if tn in used:
            continue
        return (tn, st, ovl)
    return None


def build_contig(seed, seqs, edges, target_len):
    used = set([seed])
    # path keeps (name, orient): orient +1 means forward, -1 means rc
    path = [(seed, +1)]
    # assemble sequence string; start with seed forward
    seq = seqs[seed]

    # extend to the right
    cur = seed
    cur_or = +1
    while len(seq) < target_len:
        nxt = choose_next(cur, used, edges)
        if nxt is None:
            break
        tn, st, ovl = nxt
        # orientation logic: if current is +, st='+' means tn forward appended; st='-' means tn rc
        t_or = (
            +1 if ((cur_or == +1 and st == "+") or (cur_or == -1 and st == "-")) else -1
        )
        t_seq = seqs[tn] if t_or == +1 else rc(seqs[tn])
        if ovl >= len(t_seq) or ovl >= len(seq):
            used.add(tn)
            cur = tn
            cur_or = t_or
            continue
        seq = seq + t_seq[ovl:]
        used.add(tn)
        path.append((tn, t_or))
        cur = tn
        cur_or = t_or

    # (Optional) try a short extension to the left as well for completeness
    # Reverse path and repeat once
    left_cur, left_or = path[0]
    # You could implement symmetric left-extension if desired

    return seq, [x[0] for x in path], used


def main():
    ap = argparse.ArgumentParser(
        description="Greedy ~100 kb constructor from PAF overlaps (HiFi)."
    )
    ap.add_argument("--reads", required=True, help="HiFi reads (FASTA/FASTQ, .gz ok)")
    ap.add_argument(
        "--paf", required=True, help="PAF from minimap2 all-vs-all overlaps"
    )
    ap.add_argument(
        "--min-id", type=float, default=0.985, help="min identity (fraction) [0.985]"
    )
    ap.add_argument(
        "--min-ovl", type=int, default=2000, help="min overlap length [2000]"
    )
    ap.add_argument(
        "--target", type=int, default=100000, help="target contig length [100000]"
    )
    ap.add_argument(
        "--max-contigs", type=int, default=5, help="how many contigs to attempt [5]"
    )
    ap.add_argument("--out", required=True, help="output FASTA path")
    args = ap.parse_args()

    print(f"[INFO] Loading reads from {args.reads} ...", file=sys.stderr)
    seqs = read_fx(args.reads)
    print(f"[INFO] Loaded {len(seqs)} reads.", file=sys.stderr)

    print(f"[INFO] Reading overlaps from {args.paf} ...", file=sys.stderr)
    edges = parse_paf(args.paf, args.min_id, args.min_ovl)
    print(f"[INFO] Graph nodes with edges: {len(edges)}", file=sys.stderr)

    # seed order: longest reads first
    seeds = sorted(seqs.keys(), key=lambda k: len(seqs[k]), reverse=True)

    written = 0
    used_global = set()
    with open(args.out, "w") as out:
        for s in seeds:
            if s in used_global:
                continue
            contig, used_path, used_local = build_contig(s, seqs, edges, args.target)
            if len(contig) >= min(
                args.target, int(args.target * 0.7)
            ):  # accept ≥70% target
                written += 1
                used_global |= set(used_path)
                hdr = f">unitig{written}_len{len(contig)}_seed={s}_nreads={len(used_path)}"
                out.write(hdr + "\n")
                for i in range(0, len(contig), 80):
                    out.write(contig[i : i + 80] + "\n")
                print(f"[INFO] wrote {hdr}", file=sys.stderr)
            if written >= args.max_contigs:
                break

    if written == 0:
        print(
            "[WARN] No contig reached sufficient length. Consider relaxing --min-id or --min-ovl, or recruit reads.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
