#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
polap-py-syncmer-connectivity-select-mt.py  v1.1.0

Fast, multicore organelle read selection via syncmer-connectivity PPR.

Addressing:
  --config-path FILE.yaml
  --config-dir DIR --preset NAME         # resolves DIR/NAME.yaml (default: ~/.polap/profiles)
  (Flat YAML; CLI overrides win.)

Label & Anchors (strict):
  --label mt|pt
  --mt-anchors FILE   (required if label=mt and WGS inputs)
  --pt-anchors FILE   (required if label=pt and WGS inputs)
  (No generic --anchors.)

Inputs:
  --sm-dump TSV     read_id \t hex64,hex64,...   (unique per read)
  --reads  FASTQ    (fallback; computes syncmers with -k/-s and optional HPC)

Core knobs:
  -k, -s, --hpc
  --max-occ, --min-shared, --jaccard-min, --edge-norm, --topk-nei, --steps
  --ppr-alpha, --ppr-iter, --seed-fraction
  --x-prior, --x-tsv, --nuc-cut-log10, --x-slope
  --edge-mode external|memory
  --workers, --sort-parallel, --sort-mem, --tmpdir
  --external-threshold (forces external mode on large N)
  -o, --out (out_prefix; required)
  -v/--verbose (repeatable)

Auto-tune:
  • Detects CPU & RAM to pick workers/sort-parallel/sort-mem/tmpdir.

Sharding (fixed):
  • Shard by the **last 8 hex digits** (low 32 bits).
  • Keep shard files **open** across the pass.

Outputs (labeled):
  <out>.<label>.ids
  <out>.<label>.nuclear.ids
  <out>.<label>.scores.tsv
  <out>.<label>.graph.stats.txt
"""

from __future__ import annotations
import sys, os, math, argparse, tempfile, subprocess, csv, shutil, json
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor
from typing import Dict, Any, List, Tuple, Optional

# ---------------- verbosity ----------------
VERBOSE = 0


def v(level: int, *a):
    if VERBOSE >= level:
        print("[ppr]", *a, file=sys.stderr)


# ---------------- flat YAML loader ----------------
def _load_yaml(path: str) -> Dict[str, Any]:
    if not path or not os.path.exists(path):
        return {}
    try:
        import yaml  # type: ignore

        with open(path, "r", encoding="utf-8") as fh:
            d = yaml.safe_load(fh) or {}
        if not isinstance(d, dict):
            raise ValueError("YAML root is not a mapping")
        return dict(d)
    except ImportError:
        d: Dict[str, Any] = {}
        with open(path, "r", encoding="utf-8") as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln or ln.startswith("#") or ":" not in ln:
                    continue
                k, v = ln.split(":", 1)
                k = k.strip()
                v = v.strip().strip("'").strip('"')
                vl = v.lower()
                if vl in ("true", "false"):
                    d[k] = vl == "true"
                else:
                    try:
                        if "." in v:
                            d[k] = float(v)
                        else:
                            d[k] = int(v)
                    except ValueError:
                        d[k] = v
        return d


def _resolve_cfg_path(
    config_dir: Optional[str], preset: Optional[str], config_path: Optional[str]
) -> Optional[str]:
    if config_path:
        return os.path.expanduser(config_path)
    if config_dir and preset:
        return os.path.join(os.path.expanduser(config_dir), f"{preset}.yaml")
    return None


# ---------------- small utils ----------------
def eprint(*a, **kw):
    print(*a, file=sys.stderr, **kw)


def open_guess(fn):
    import gzip

    return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn, "rt")


def norm_id(s: str) -> str:
    return (s or "").split()[0]


def bytes_to_human(n):
    if n is None:
        return "unknown"
    for u in ["B", "KB", "MB", "GB", "TB"]:
        if n < 1024:
            return f"{n:.1f} {u}"
        n /= 1024
    return f"{n:.1f} PB"


def _read_total_ram():
    try:
        import psutil  # type: ignore

        return int(psutil.virtual_memory().total)
    except Exception:
        pass
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    kB = int(line.split()[1])
                    return kB * 1024
    except Exception:
        return None
    return None


def _pick_tmpdir(user_tmpdir=None, need_bytes=0):
    if user_tmpdir:
        return user_tmpdir
    shm = "/dev/shm"
    try:
        if os.path.isdir(shm):
            total, used, free = shutil.disk_usage(shm)
            if free >= need_bytes:
                return shm
    except Exception:
        pass
    return tempfile.gettempdir()


def _pick_sort_bin():
    return shutil.which("gsort") or shutil.which("sort") or "sort"


SORT_BIN = _pick_sort_bin()


# ---------------- FASTQ + syncmers (fallback path) ----------------
def fastq_iter(path: str):
    with open_guess(path) as fh:
        while True:
            n = fh.readline()
            if not n:
                break
            s = fh.readline()
            fh.readline()
            q = fh.readline()
            if not q:
                break
            yield norm_id(n[1:].strip()), s.strip()


def nt4(c: str) -> int:
    if c in "Aa":
        return 0
    if c in "Cc":
        return 1
    if c in "Gg":
        return 2
    if c in "Tt":
        return 3
    return -1


def rc2b(x: int, l: int) -> int:
    mask = (1 << (2 * l)) - 1 if l < 32 else (1 << 64) - 1
    x = (~x) & mask
    y = 0
    for _ in range(l):
        y = (y << 2) | (x & 3)
        x >>= 2
    return y


def hpc_compress(seq: str) -> str:
    out = []
    prev = ""
    for c in seq.upper():
        if c not in "ACGT":
            prev = ""
            continue
        if c != prev:
            out.append(c)
            prev = c
    return "".join(out)


def closed_syncmers(seq: str, k: int, s: int):
    n = len(seq)
    if n < k:
        return []
    ns = n - s + 1
    w = k - s + 1
    mask = (1 << (2 * s)) - 1 if s < 32 else (1 << 64) - 1
    f = r = l = 0
    sh = [0] * ns
    nw = 0
    for c in seq:
        b = nt4(c)
        if b < 0:
            l = f = r = 0
            continue
        f = ((f << 2) | b) & mask
        r = (rc2b(b, 1) << (2 * (s - 1))) | (r >> 2)
        if l < s:
            l += 1
            if l < s:
                continue
        can = f if f < r else r
        if nw < ns:
            sh[nw] = can
            nw += 1
    if nw > ns:
        nw = ns
    out = []
    for i in range(0, n - k + 1):
        mn = (1 << 64) - 1
        arg = -1
        for j in range(w):
            v = sh[i + j]
            if v < mn:
                mn = v
                arg = j
        if arg == 0 or arg == w - 1:
            out.append(mn)
    return out


# ---------------- anchors & X ----------------
def load_id_set(path: Optional[str]) -> set[str]:
    S = set()
    if not path:
        return S
    with open(path) as fh:
        for line in fh:
            t = line.strip().split()
            if t:
                S.add(t[0])
    return S


def load_X_quickview(tsv: Optional[str]) -> Dict[str, float]:
    X = {}
    if not tsv:
        return X
    with open(tsv, newline="") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            try:
                rid = norm_id(row["read_id"])
                med = float(row["median"])
                if med > 0:
                    X[rid] = math.log10(med)
            except:
                pass
    return X


# ---------------- external helpers ----------------
def _sort_uniq_count(in_path, out_path, parallel=None, mem=None):
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    cmd = [SORT_BIN, in_path]
    if parallel or mem:
        cmd = (
            [SORT_BIN]
            + ([f"--parallel={parallel}"] if parallel else [])
            + ([f"-S{mem}"] if mem else [])
            + [in_path]
        )
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, env=env, text=True)
    p2 = subprocess.Popen(
        ["uniq", "-c"], stdin=p1.stdout, stdout=open(out_path, "w"), env=env, text=True
    )
    p1.stdout.close()
    p2.communicate()


def _sort_key1(in_path, out_path, parallel=None, mem=None):
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    cmd = [SORT_BIN, "-k1,1", in_path]
    if parallel or mem:
        cmd = (
            [SORT_BIN, "-k1,1"]
            + ([f"--parallel={parallel}"] if parallel else [])
            + ([f"-S{mem}"] if mem else [])
            + [in_path]
        )
    subprocess.check_call(cmd, stdout=open(out_path, "w"), env=env, text=True)


def _sort_key12_numeric(in_path, out_path, parallel=None, mem=None):
    env = os.environ.copy()
    env["LC_ALL"] = "C"
    cmd = [SORT_BIN, "-k1,1n", "-k2,2n", in_path]
    if parallel or mem:
        cmd = (
            [SORT_BIN, "-k1,1n", "-k2,2n"]
            + ([f"--parallel={parallel}"] if parallel else [])
            + ([f"-S{mem}"] if mem else [])
            + [in_path]
        )
    subprocess.check_call(cmd, stdout=open(out_path, "w"), env=env, text=True)


# ---------------- sm-dump loader (streaming/external) ----------------
def load_sm_dump_external_stream(
    path, max_occ, tmpdir=None, sort_parallel=None, sort_mem=None
):
    """Return (read_index, idx2read, card, sm_read_filtered_sorted_path) with max_occ filter."""
    v(1, "[load] pass1: assign indices & spill all hashes and sm_read …")
    read_index = {}
    idx2read = []
    card = []

    all_fd, all_path = tempfile.mkstemp(prefix="sm_all_", dir=tmpdir, text=True)
    os.close(all_fd)
    smr_fd, smr_path = tempfile.mkstemp(prefix="sm_read_", dir=tmpdir, text=True)
    os.close(smr_fd)

    with open(path) as fh, open(all_path, "w") as wh, open(smr_path, "w") as wr:
        for n, line in enumerate(fh, 1):
            if (n % 1_000_000) == 0:
                v(1, f"[load] pass1 … {n/1e6:.1f}M lines")
            line = line.rstrip("\n")
            if not line:
                continue
            rid, hashes = (line.split("\t", 1) + [""])[:2]
            rid = rid.strip()
            if not rid:
                continue
            if rid in read_index:
                continue
            vals = []
            if hashes:
                for x in hashes.split(","):
                    if x:
                        vals.append(int(x, 16))
            if len(vals) > 1:
                vals = sorted(set(vals))
            idx = len(idx2read)
            read_index[rid] = idx
            idx2read.append(rid)
            card.append(len(vals))
            for h in vals:
                wh.write(f"{h:016x}\n")
                wr.write(f"{h:016x}\t{idx}\n")

    v(1, "[load] counting hashes (sort|uniq -c) …")
    cnt_fd, cnt_path = tempfile.mkstemp(prefix="sm_cnt_", dir=tmpdir, text=True)
    os.close(cnt_fd)
    _sort_uniq_count(all_path, cnt_path, parallel=sort_parallel, mem=sort_mem)
    os.remove(all_path)

    v(1, "[load] derive allowed-hash list …")
    allow_fd, allow_path = tempfile.mkstemp(prefix="sm_allow_", dir=tmpdir, text=True)
    os.close(allow_fd)
    with open(cnt_path) as fi, open(allow_path, "w") as fo:
        for line in fi:
            parts = line.strip().split()
            if len(parts) != 2:
                continue
            occ = int(parts[0])
            hx = parts[1]
            if occ <= max_occ:
                fo.write(hx + "\n")
    os.remove(cnt_path)

    v(1, "[load] sort sm_read by hash …")
    smr_sorted_fd, smr_sorted_path = tempfile.mkstemp(
        prefix="sm_read_sorted_", dir=tmpdir, text=True
    )
    os.close(smr_sorted_fd)
    _sort_key1(smr_path, smr_sorted_path, parallel=sort_parallel, mem=sort_mem)
    os.remove(smr_path)

    v(1, "[load] filter sm_read by allowed hashes (sorted merge) …")
    filt_fd, filt_path = tempfile.mkstemp(
        prefix="sm_read_allowed_", dir=tmpdir, text=True
    )
    os.close(filt_fd)
    with open(allow_path) as fa, open(smr_sorted_path) as fs, open(
        filt_path, "w"
    ) as fo:
        ha = fa.readline().strip()
        hs = fs.readline().strip()
        while hs:
            if not ha:
                break
            hx = hs.split("\t", 1)[0]
            if hx == ha:
                fo.write(hs + "\n")
                hs = fs.readline().strip()
            elif hx < ha:
                hs = fs.readline().strip()
            else:
                ha = fa.readline().strip()
    os.remove(allow_path)
    os.remove(smr_sorted_path)
    return read_index, idx2read, card, filt_path


# ---------------- sharded edge build ----------------
def _shard_index(hex_hash: str, n_shards: int) -> int:
    tail = hex_hash[-8:] if len(hex_hash) >= 8 else hex_hash
    try:
        x = int(tail, 16)
    except Exception:
        x = 0
    return x % max(n_shards, 1)


def _pairs_aggregate_one(
    shard_path: str, out_counts_path: str, sort_parallel=None, sort_mem=None
):
    pair_fd, pair_path = tempfile.mkstemp(
        prefix="pairs_", dir=os.path.dirname(out_counts_path), text=True
    )
    os.close(pair_fd)
    with open(shard_path) as f, open(pair_path, "w") as wp:
        w = wp.write
        cur_hash = None
        group = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            hx, idx_s = line.split("\t")
            idx = int(idx_s)
            if cur_hash is None:
                cur_hash = hx
                group = [idx]
            elif hx == cur_hash:
                group.append(idx)
            else:
                if len(group) >= 2:
                    group.sort()
                    L = len(group)
                    for i in range(L):
                        u = group[i]
                        for j in range(i + 1, L):
                            v = group[j]
                            w(f"{u}\t{v}\n")
                cur_hash = hx
                group = [idx]
        if group and len(group) >= 2:
            group.sort()
            L = len(group)
            for i in range(L):
                u = group[i]
                for j in range(i + 1, L):
                    v = group[j]
                    w(f"{u}\t{v}\n")
    counts_tmp_fd, counts_tmp = tempfile.mkstemp(
        prefix="pairs_cnt_", dir=os.path.dirname(out_counts_path), text=True
    )
    os.close(counts_tmp_fd)
    _sort_uniq_count(pair_path, counts_tmp, parallel=sort_parallel, mem=sort_mem)
    os.remove(pair_path)
    with open(counts_tmp) as fi, open(out_counts_path, "w") as fo:
        for line in fi:
            parts = line.strip().split()
            if len(parts) != 3:
                continue
            c, u, v = parts[0], parts[1], parts[2]
            fo.write(f"{u}\t{v}\t{c}\n")
    os.remove(counts_tmp)
    return out_counts_path


def build_edges_external_stream_sharded(
    sm_read_filtered_path,
    workers,
    min_shared,
    jacc_min,
    edge_norm,
    topk_nei,
    card,
    tmpdir=None,
    sort_parallel=None,
    sort_mem=None,
):
    if not workers or workers < 1:
        workers = max(os.cpu_count() or 1, 1)

    # 1) shard by hash (keep files open)
    sh_paths = []
    outs = []
    for i in range(workers):
        fd, p = tempfile.mkstemp(prefix=f"sh_{i:02d}_", dir=tmpdir, text=True)
        os.close(fd)
        sh_paths.append(p)
        outs.append(open(p, "w"))
    try:
        with open(sm_read_filtered_path) as f:
            for line in f:
                hx = line.split("\t", 1)[0]
                s = _shard_index(hx, workers)
                outs[s].write(line)
    finally:
        for w in outs:
            w.close()

    # 2) per-shard aggregate in parallel
    shard_outs = []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = []
        for p in sh_paths:
            fd, outp = tempfile.mkstemp(prefix="cnt_", dir=tmpdir, text=True)
            os.close(fd)
            futs.append(
                (
                    ex.submit(_pairs_aggregate_one, p, outp, sort_parallel, sort_mem),
                    p,
                    outp,
                )
            )
        for fut, p, outp in futs:
            _ = fut.result()
            shard_outs.append(outp)
            try:
                os.remove(p)
            except:
                pass

    # 3) merged/sorted
    merged_fd, merged_path = tempfile.mkstemp(prefix="cnt_all_", dir=tmpdir, text=True)
    os.close(merged_fd)
    with open(merged_path, "w") as w:
        for p in shard_outs:
            with open(p) as r:
                for line in r:
                    w.write(line)
            try:
                os.remove(p)
            except:
                pass

    sorted_fd, sorted_path = tempfile.mkstemp(
        prefix="cnt_all_sorted_", dir=tmpdir, text=True
    )
    os.close(sorted_fd)
    _sort_key12_numeric(merged_path, sorted_path, parallel=sort_parallel, mem=sort_mem)
    os.remove(merged_path)

    # 4) fold + build adjacency
    adj = {}
    with open(sorted_path) as f:
        prev_u = prev_v = None
        total = 0
        for line in f:
            u_s, v_s, c_s = line.strip().split("\t")
            u = int(u_s)
            v = int(v_s)
            c = int(c_s)
            if prev_u is None:
                prev_u, prev_v, total = u, v, c
                continue
            if u == prev_u and v == prev_v:
                total += c
                continue
            _flush_pair(
                adj, prev_u, prev_v, total, min_shared, jacc_min, edge_norm, card
            )
            prev_u, prev_v, total = u, v, c
        if prev_u is not None:
            _flush_pair(
                adj, prev_u, prev_v, total, min_shared, jacc_min, edge_norm, card
            )
    os.remove(sorted_path)

    if topk_nei and topk_nei > 0:
        for u, nbrs in list(adj.items()):
            if len(nbrs) > topk_nei:
                nbrs.sort(key=lambda x: x[1], reverse=True)
                adj[u] = nbrs[:topk_nei]
    return adj


def _flush_pair(adj, u, v, inter, min_shared, jacc_min, edge_norm, card):
    if inter < min_shared:
        return
    union = card[u] + card[v] - inter
    if union <= 0:
        return
    jacc = inter / union
    if jacc < jacc_min:
        return
    w = inter
    if edge_norm:
        w = inter / max(1.0, math.sqrt(card[u] * card[v]))
    adj.setdefault(u, []).append((v, w))
    adj.setdefault(v, []).append((u, w))


# ---------------- PPR & threshold ----------------
def bfs_expand(adj, seeds, steps):
    cur = set(seeds)
    frontier = set(seeds)
    for _ in range(steps):
        nxt = set()
        for u in frontier:
            for v, _ in adj.get(u, []):
                if v not in cur:
                    nxt.add(v)
        if not nxt:
            break
        cur |= nxt
        frontier = nxt
    return cur


def personalized_pagerank(adj, seeds, alpha=0.85, iters=30, seed_fraction=1.0):
    nodes = sorted(adj.keys())
    ix = {u: i for i, u in enumerate(nodes)}
    N = len(nodes)
    if N == 0:
        return {}
    nbr = [adj.get(u, [])[:] for u in nodes]
    for i in range(N):
        tot = sum(w for _, w in nbr[i]) or 1.0
        nbr[i] = [(v, w / tot) for (v, w) in nbr[i]]
    s = [0.0] * N
    if seeds:
        mass = seed_fraction / len(seeds)
        for u in seeds:
            if u in ix:
                s[ix[u]] += mass
    p = s[:]
    for _ in range(iters):
        new = [0.0] * N
        for i in range(N):
            if nbr[i]:
                pi = p[i]
                for v, w in nbr[i]:
                    j = ix.get(v, None)
                    if j is not None:
                        new[j] += alpha * pi * w
        for j in range(N):
            new[j] += (1 - alpha) * s[j]
        p = new
    return {nodes[i]: p[i] for i in range(N)}


def apply_x_prior(scores, idx2read, Xmap, nuc_cut_log10, x_slope):
    if nuc_cut_log10 is None:
        return scores

    def sigmoid(z):
        return 1.0 / (1.0 + math.exp(-z))

    slope = max(float(x_slope), 1e-6)
    out = {}
    for u, sc in scores.items():
        rid = idx2read[u]
        x = Xmap.get(rid, nuc_cut_log10)
        out[u] = sc * sigmoid((x - nuc_cut_log10) / slope)
    return out


def _auto_threshold(values):
    try:
        import numpy as np

        x = np.array(values, dtype=float)
        if x.size == 0:
            return 0.0
        hist, edges = np.histogram(x, bins=64)
        if hist.sum() == 0:
            return float(np.percentile(x, 80))
        w0 = np.cumsum(hist) / hist.sum()
        centers = 0.5 * (edges[:-1] + edges[1:])
        mu = np.cumsum(hist * centers) / hist.sum()
        mu_t = mu[-1] if mu.size else 0.0
        w1 = 1.0 - w0
        between = (mu_t * w0 - mu) ** 2 / (w0 * w1 + 1e-12)
        k = int(np.nanargmax(between)) if between.size else int(len(x) * 0.2)
        return float(edges[min(k + 1, len(edges) - 1)])
    except Exception:
        v = sorted(values)
        if not v:
            return 0.0
        pos = int(0.8 * (len(v) - 1))
        return float(v[pos])


# ---------------- auto-tune ----------------
def _auto_tune(args):
    cpu = max(os.cpu_count() or 1, 1)
    ram = _read_total_ram() or 8 * 1024**3
    need_tmp = int(0.2 * ram)
    if not args.tmpdir:
        args.tmpdir = _pick_tmpdir(None, need_tmp)
    if not args.workers or args.workers <= 0:
        if cpu <= 4:
            args.workers = 2
        elif cpu <= 8:
            args.workers = 4
        elif cpu <= 16:
            args.workers = 6
        elif cpu <= 32:
            args.workers = 12
        else:
            args.workers = min(16, cpu // 2)
    if not args.sort_parallel or args.sort_parallel <= 0:
        sp = max(1, cpu // max(1, args.workers))
        if cpu >= 8 and sp < 2:
            sp = 2
        args.sort_parallel = sp
    if not args.sort_mem:
        reserve = 0.15
        mem_for_sort = int(ram * (1.0 - reserve) * 0.6)
        max_threads = max(1, args.workers * args.sort_parallel)
        per = max(512 * 1024**2, mem_for_sort // max_threads)
        per = min(per, 6 * 1024**3)

        def _S(n):
            if n >= 1024**3:
                return f"{n//(1024**3)}G"
            if n >= 1024**2:
                return f"{n//(1024**2)}M"
            return f"{max(1,n//1024)}K"

        args.sort_mem = _S(per)
    v(1, f"[auto] cpu={cpu} ram={bytes_to_human(ram)} tmpdir={args.tmpdir}")
    v(
        1,
        f"[auto] workers={args.workers} sort-parallel={args.sort_parallel} sort-mem={args.sort_mem}",
    )


# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser(
        description="Syncmer-connectivity PPR (organelle selection) — auto-tuned multicore external graph build."
    )
    # addressing
    ap.add_argument("--config-path")
    ap.add_argument("--config-dir")
    ap.add_argument("--preset")
    # verbosity
    ap.add_argument("-v", "--verbose", action="count", default=0)
    # label + anchors (strict)
    ap.add_argument("--label", default="mt", choices=["mt", "pt"])
    ap.add_argument("--mt-anchors", help="anchor ID list for mt")
    ap.add_argument("--pt-anchors", help="anchor ID list for pt")
    # inputs
    ap.add_argument(
        "--sm-dump", help="TSV: read_id\\thex64,hex64,... (unique per read)"
    )
    ap.add_argument("--reads", help="FASTQ(.gz) if no sm-dump")
    # syncmer params for FASTQ fallback
    ap.add_argument("-k", type=int)
    ap.add_argument("-s", type=int)
    ap.add_argument("--hpc", action="store_true")
    # graph knobs
    ap.add_argument("--max-occ", type=int)
    ap.add_argument("--min-shared", type=int)
    ap.add_argument("--jaccard-min", type=float)
    ap.add_argument("--edge-norm", action="store_true")
    ap.add_argument("--topk-nei", type=int)
    ap.add_argument("--steps", type=int)
    # ppr
    ap.add_argument("--ppr-alpha", type=float)
    ap.add_argument("--ppr-iter", type=int)
    ap.add_argument("--seed-fraction", type=float)
    # x-prior
    ap.add_argument("--x-prior", action="store_true")
    ap.add_argument("--x-tsv")
    ap.add_argument("--nuc-cut-log10", type=float)
    ap.add_argument("--x-slope", type=float)
    ap.add_argument(
        "--score-th", type=float, help="threshold for final PPR score (<=0 -> auto)"
    )
    # engine
    ap.add_argument("--edge-mode", choices=["external", "memory"])
    ap.add_argument("--workers", type=int)
    ap.add_argument("--tmpdir")
    ap.add_argument("--external-threshold", type=int)
    ap.add_argument("--sort-parallel", type=int)
    ap.add_argument("--sort-mem")
    # output
    ap.add_argument("-o", "--out", required=True, help="out_prefix for labeled outputs")
    ap.add_argument("--version", action="store_true")

    args = ap.parse_args()

    global VERBOSE
    VERBOSE = int(args.verbose or 0)

    if args.version:
        print("1.1.0")
        return 0

    # YAML load + merge (CLI wins)
    cfg_path = _resolve_cfg_path(args.config_dir, args.preset, args.config_path)
    cfg = _load_yaml(cfg_path) if cfg_path else {}
    v(1, f"[cfg] config={cfg_path or '(none)'}")

    def val(key, default=None):
        # map flat YAML key to attribute name (replace '_' with '-' on CLI; we read CLI first)
        cli_name = key.replace("_", "-")
        # prefer exact attribute if present in argparse
        attr = cli_name.replace("-", "_")
        cli = getattr(args, attr, None)
        if cli is not None:
            return cli
        return cfg.get(key, default)

    # resolve label + anchors
    label = (val("label", "mt") or "mt").strip()
    if label not in ("mt", "pt"):
        eprint("[error] --label must be mt or pt")
        return 2
    mt_anch = val("mt_anchors")
    pt_anch = val("pt_anchors")
    if label == "mt":
        anchors = load_id_set(mt_anch)
        if not anchors:
            eprint("[error] no anchors for label=mt (mt_anchors or --mt-anchors)")
            return 2
    else:
        anchors = load_id_set(pt_anch)
        if not anchors:
            eprint("[error] no anchors for label=pt (pt_anchors or --pt-anchors)")
            return 2

    # inputs: sm-dump or reads
    sm_dump = val("sm_dump")
    reads = val("reads")
    if not sm_dump and not reads:
        eprint("[error] need --sm-dump or --reads")
        return 2

    # syncmer params for FASTQ fallback
    preset = (val("preset", "hifi") or "hifi").lower()
    k = int(val("k", 41 if preset == "ont" else 121))
    s = int(val("s", 21 if preset == "ont" else 27))
    hpc = bool(val("hpc", preset == "ont"))

    # graph knobs (defaults)
    max_occ = int(val("max_occ", 200))
    min_shared = int(val("min_shared", 5 if preset == "ont" else 4))
    jacc_min = float(val("jaccard_min", 0.0075 if preset == "ont" else 0.015))
    edge_norm = bool(val("edge_norm", True))
    topk_nei = int(val("topk_nei", 40 if preset == "ont" else 50))
    steps = int(val("steps", 1 if preset == "ont" else 2))

    # ppr knobs
    ppr_alpha = float(val("ppr_alpha", 0.85))
    ppr_iter = int(val("ppr_iter", 30))
    seed_fraction = float(val("seed_fraction", 1.0))

    # x-prior knobs
    x_prior = bool(val("x_prior", False))
    x_tsv = val("x_tsv")
    nuc_cut_log10 = val("nuc_cut_log10", None)
    nuc_cut_log10 = float(nuc_cut_log10) if nuc_cut_log10 is not None else None
    x_slope = float(val("x_slope", 0.20))

    # engine knobs
    edge_mode = (val("edge_mode", "external") or "external").lower()
    workers = val("workers", None)
    workers = int(workers) if workers not in (None, "") else None
    tmpdir = val("tmpdir", None)
    external_threshold = int(val("external_threshold", 500000))
    sort_parallel = val("sort_parallel", None)
    sort_parallel = int(sort_parallel) if sort_parallel not in (None, "") else None
    sort_mem = val("sort_mem", None)

    out_prefix = val("out", None) or val("out_prefix", None)
    if not out_prefix:
        eprint("[error] need -o/--out (out_prefix)")
        return 2
    outdir = os.path.dirname(out_prefix) or "."
    os.makedirs(outdir, exist_ok=True)

    # auto-tune
    class Dummy:
        pass

    a = Dummy()
    a.workers = workers
    a.sort_parallel = sort_parallel
    a.sort_mem = sort_mem
    a.tmpdir = tmpdir
    _auto_tune(a)
    workers, sort_parallel, sort_mem, tmpdir = (
        a.workers,
        a.sort_parallel,
        a.sort_mem,
        a.tmpdir,
    )

    # ingest
    if sm_dump:
        v(1, "[step] loading precomputed syncmers (external stream) …")
        read_index, idx2read, card, sm_read_filtered_path = (
            load_sm_dump_external_stream(
                sm_dump,
                max_occ,
                tmpdir=tmpdir,
                sort_parallel=sort_parallel,
                sort_mem=sort_mem,
            )
        )
    else:
        v(1, "[step] computing syncmers from FASTQ (fallback; slower) …")
        sm_counts = Counter()
        per_read_sms = []
        idx2read = []
        read_index = {}
        for rid, seq in fastq_iter(reads):
            if hpc:
                seq = hpc_compress(seq)
            sms = closed_syncmers(seq, k, s)
            sms = sorted(set(sms)) if sms else []
            read_index[rid] = len(idx2read)
            idx2read.append(rid)
            per_read_sms.append(sms)
            for h in sms:
                sm_counts[h] += 1
        card = [0] * len(idx2read)
        fd, path = tempfile.mkstemp(prefix="sm_read_allowed_", dir=tmpdir, text=True)
        os.close(fd)
        with open(path, "w") as fo:
            for i, sms in enumerate(per_read_sms):
                card[i] = len(sms)
                if not sms:
                    continue
                keep = [h for h in sms if sm_counts[h] <= max_occ]
                keep = sorted(set(keep))
                for h in keep:
                    fo.write(f"{h:016x}\t{i}\n")
        sorted_fd, sm_read_filtered_path = tempfile.mkstemp(
            prefix="sm_read_filt_sorted_", dir=tmpdir, text=True
        )
        os.close(sorted_fd)
        _sort_key1(path, sm_read_filtered_path, parallel=sort_parallel, mem=sort_mem)
        os.remove(path)

    N = len(idx2read)
    v(1, f"[info] reads_in_dump={N}")
    # graph build mode
    mode = edge_mode
    if mode == "memory" and N > external_threshold:
        mode = "external"
    v(1, f"[step] building graph mode={mode} …")

    if mode == "external":
        adj = build_edges_external_stream_sharded(
            sm_read_filtered_path,
            workers,
            min_shared,
            jacc_min,
            edge_norm,
            topk_nei,
            card,
            tmpdir=tmpdir,
            sort_parallel=sort_parallel,
            sort_mem=sort_mem,
        )
    else:
        # in-memory fallback
        sm2reads = defaultdict(list)
        with open(sm_read_filtered_path) as f:
            for line in f:
                hx, idx_s = line.strip().split("\t")
                sm2reads[hx].append(int(idx_s))
        adj = build_edges_memory(
            sm2reads, min_shared, jacc_min, edge_norm, topk_nei, card
        )
    try:
        os.remove(sm_read_filtered_path)
    except:
        pass

    # seeds
    seed_idx = [read_index[r] for r in anchors if r in read_index]
    if not seed_idx:
        eprint("[error] none of the anchors found in dump")
        return 2

    # BFS -> subgraph
    v(1, "[step] BFS candidate expansion …")
    cand = bfs_expand(adj, seed_idx, steps)
    sub_adj = {
        u: [(v, w) for (v, w) in adj.get(u, []) if v in cand] for u in cand if u in adj
    }
    v(1, f"[info] cand_nodes={len(sub_adj)}")

    # PPR
    v(1, "[step] PPR …")
    scores = personalized_pagerank(
        sub_adj, seed_idx, alpha=ppr_alpha, iters=ppr_iter, seed_fraction=seed_fraction
    )

    # X-prior
    if x_prior and nuc_cut_log10 is not None and x_tsv:
        v(1, "[step] applying X-prior …")
        Xmap = load_X_quickview(x_tsv)
        scores = apply_x_prior(scores, idx2read, Xmap, nuc_cut_log10, x_slope or 0.20)

    # threshold
    vals = list(scores.values())
    th = float(val("score_th", 0.0) or 0.0)
    if th <= 0.0:
        th = _auto_threshold(vals)
    v(1, f"[info] threshold={th:.6g}")

    # outputs
    os.makedirs(os.path.dirname(out_prefix) or ".", exist_ok=True)
    ids_path = f"{out_prefix}.{label}.ids"
    nuc_path = f"{out_prefix}.{label}.nuclear.ids"
    sc_path = f"{out_prefix}.{label}.scores.tsv"
    st_path = f"{out_prefix}.{label}.graph.stats.txt"

    with open(ids_path, "w") as w_ids, open(nuc_path, "w") as w_nuc, open(
        sc_path, "w"
    ) as w_sc:
        w_sc.write("read_id\tdeg\tshared_sum\tfinal_score\n")
        for u in sorted(sub_adj.keys()):
            rid = idx2read[u]
            deg = len(sub_adj[u])
            shared_sum = sum(w for _, w in sub_adj[u])
            s = scores.get(u, 0.0)
            w_sc.write(f"{rid}\t{deg}\t{shared_sum:.3f}\t{s:.6g}\n")
            (w_ids if s >= th else w_nuc).write(rid + "\n")

    with open(st_path, "w") as w:
        w.write(f"nodes\t{len(adj)}\n")
        w.write(f"cand_nodes\t{len(sub_adj)}\n")
        w.write(f"anchors\t{len(seed_idx)}\n")
        w.write(f"max_occ\t{max_occ}\n")
        w.write(f"min_shared\t{min_shared}\n")
        w.write(f"jaccard_min\t{jacc_min}\n")
        w.write(f"edge_norm\t{int(bool(edge_norm))}\n")
        w.write(f"topk_nei\t{topk_nei}\n")
        w.write(f"steps\t{steps}\n")
        w.write(f"ppr_alpha\t{ppr_alpha}\n")
        w.write(f"ppr_iter\t{ppr_iter}\n")
        w.write(
            f"x_prior\t{int(bool(x_prior and nuc_cut_log10 is not None and x_tsv))}\n"
        )
        if x_prior and nuc_cut_log10 is not None:
            w.write(f"nuc_cut_log10\t{nuc_cut_log10}\n")
            w.write(f"x_slope\t{x_slope}\n")
        w.write(f"threshold\t{th}\n")
        w.write(f"label\t{label}\n")

    v(1, f"[done] {ids_path}")
    v(1, f"[done] {nuc_path}")
    v(1, f"[done] {sc_path}")
    v(1, f"[done] {st_path}")
    return 0


# ---------------- in-memory fallback ----------------
def build_edges_memory(sm2reads, min_shared, jacc_min, edge_norm, topk_nei, card):
    inter = defaultdict(int)
    for lst in sm2reads.values():
        lst = sorted(lst)
        L = len(lst)
        for i in range(L):
            u = lst[i]
            for j in range(i + 1, L):
                v = lst[j]
                inter[(u, v)] += 1
    adj = {}
    for (u, v), cnt in inter.items():
        if cnt < min_shared:
            continue
        union = card[u] + card[v] - cnt
        if union <= 0:
            continue
        jacc = cnt / union
        if jacc < jacc_min:
            continue
        w = cnt
        if edge_norm:
            w = cnt / max(1.0, math.sqrt(card[u] * card[v]))
        adj.setdefault(u, []).append((v, w))
        adj.setdefault(v, []).append((u, w))
    if topk_nei and topk_nei > 0:
        for u, nbrs in list(adj.items()):
            if len(nbrs) > topk_nei:
                nbrs.sort(key=lambda x: x[1], reverse=True)
                adj[u] = nbrs[:topk_nei]
    return adj


if __name__ == "__main__":
    sys.exit(main())
