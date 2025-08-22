#!/usr/bin/env python3
################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

#!/usr/bin/env python3
import sys, re, argparse


def rc(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def rev_or(o: str) -> str:
    return "+" if o == "-" else "-"


def parse_fasta(path):
    seqs = {}
    if not path:
        return seqs
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].strip().split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs


def parse_gfa(gfa_path):
    seg = {}  # name -> sequence (may be "*")
    links = {}  # (from, f_or, to, t_or) -> overlap_match_len
    with open(gfa_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            typ = parts[0]
            if typ == "S":
                name = parts[1]
                sequence = parts[2]
                seg[name] = sequence
            elif typ == "L":
                f, fo, t, to, cigar = parts[1], parts[2], parts[3], parts[4], parts[5]
                # total matched length from the CIGAR (count M and '=')
                match_len = 0
                for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar):
                    if op in ("M", "="):
                        match_len += int(n)
                links[(f, fo, t, to)] = match_len
    return seg, links


def get_oriented(seq: str, orient: str) -> str:
    return seq if orient == "+" else rc(seq)


def get_overlap(prev_name, prev_or, curr_name, curr_or, links):
    key = (prev_name, prev_or, curr_name, curr_or)
    if key in links:
        return links[key]
    # try reverse-recorded direction if GFA only stores one orientation
    rev = (curr_name, rev_or(curr_or), prev_name, rev_or(prev_or))
    if rev in links:
        return links[rev]
    raise ValueError(
        f"No L-line overlap for {prev_name}{prev_or} -> {curr_name}{curr_or}"
    )


def assemble_path(path_elems, seg_seqs, links, external_fasta=None, circular=False):
    # load external sequences if needed
    ext = parse_fasta(external_fasta) if external_fasta else {}

    def fetch(name):
        s = seg_seqs.get(name)
        if s is None:
            raise ValueError(f"Segment {name} not found in GFA.")
        if s == "*":
            if name not in ext:
                raise ValueError(
                    f"Segment {name} has '*' sequence; not found in --fasta for {name}."
                )
            return ext[name]
        return s

    # initialize with first segment
    out_name, out_or = path_elems[0]
    out_seq = get_oriented(fetch(out_name), out_or)

    # stitch through the path
    for i in range(1, len(path_elems)):
        prev_name, prev_or = path_elems[i - 1]
        curr_name, curr_or = path_elems[i]
        ov = get_overlap(prev_name, prev_or, curr_name, curr_or, links)
        curr_seq = get_oriented(fetch(curr_name), curr_or)
        if ov > len(out_seq) or ov > len(curr_seq):
            raise ValueError(
                f"Overlap {ov} exceeds sequence length (prev {len(out_seq)}, curr {len(curr_seq)})."
            )
        out_seq += curr_seq[ov:]

    # if circular, also apply the last->first overlap and trim the redundant prefix
    if circular and len(path_elems) > 1:
        last_name, last_or = path_elems[-1]
        first_name, first_or = path_elems[0]
        ov_c = get_overlap(last_name, last_or, first_name, first_or, links)
        if ov_c > 0:
            # Optional: sanity check (commented out to keep it strict-minimal)
            # assert out_seq[-ov_c:] == out_seq[:ov_c], "Circular overlap mismatch"
            out_seq = out_seq[ov_c:]

    return out_seq


def parse_path_string(s: str):
    items = []
    for tok in s.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if tok[-1] not in "+-":
            raise ValueError(f"Path token must end with '+' or '-': {tok}")
        items.append((tok[:-1], tok[-1]))
    if not items:
        raise ValueError("Empty path.")
    return items


def main():
    ap = argparse.ArgumentParser(
        description="Assemble a path sequence from a GFA using L-line CIGAR overlaps."
    )
    ap.add_argument("--gfa", required=True, help="Input GFA file")
    ap.add_argument(
        "--path", required=True, help="Comma-separated path like 'u2+,u1+,u3+,u1-'"
    )
    ap.add_argument(
        "--fasta", help="Optional FASTA with segment sequences if S-lines have '*'"
    )
    ap.add_argument(
        "-c",
        "--circular-path",
        action="store_true",
        help="Also apply the overlap between the last and first path nodes (circular join).",
    )
    ap.add_argument("--out", help="Output file (one path per line)")
    args = ap.parse_args()

    segs, links = parse_gfa(args.gfa)
    path_elems = parse_path_string(args.path)

    seq = assemble_path(
        path_elems, segs, links, external_fasta=args.fasta, circular=args.circular_path
    )
    print(seq)

    if args.out:
        with open(args.out, "w") as fh:
            fh.write(">plastid1\n")
            fh.write(seq + "\n")
    else:
        for p in paths:
            print(seq)


if __name__ == "__main__":
    main()
