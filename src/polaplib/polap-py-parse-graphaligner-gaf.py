#!/usr/bin/env python3
import sys
import csv
import gzip


def open_maybe_gzip(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def main():
    if len(sys.argv) < 7:
        sys.stderr.write(
            "Usage: polap-py-parse-graphaligner-gaf.py "
            "<in.gaf[.gz]> <reject.txt> <retain.txt> <summary.tsv> "
            "<identity_thresh> <clip_thresh> [min_aln_len]\n"
        )
        sys.exit(2)

    # === Inputs ===
    gaf_path = sys.argv[1]  # input GAF(.gz)
    reject_path = sys.argv[2]  # output: rejected read names
    retain_path = sys.argv[3]  # output: retained read names
    summary_path = sys.argv[4]  # output: summary.tsv
    identity_thr = float(sys.argv[5])  # e.g., 0.95
    clip_thr = int(sys.argv[6])  # e.g., 100
    min_aln_len = int(sys.argv[7]) if len(sys.argv) >= 8 else 0  # e.g., 1000

    # === Open files ===
    reject_out = open(reject_path, "w")
    retain_out = open(retain_path, "w")
    summary_out = open(summary_path, "w", newline="")
    summary_writer = csv.writer(summary_out, delimiter="\t")
    summary_writer.writerow(
        ["read_id", "identity", "left_clip", "right_clip", "aln_len", "keep"]
    )

    # GAF columns (0-based):
    # 0 qname, 1 qlen, 2 qstart, 3 qend, 4 strand, 5 tname, 6 tlen, 7 tstart, 8 tend,
    # 9 matches, 10 alnBlockLen, 11 mapq, ...
    with open_maybe_gzip(gaf_path, "rt") as f:
        for line in f:
            if not line or line[0] == "#":
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue  # not a valid GAF line

            try:
                qname = fields[0]
                qlen = int(fields[1])
                qstart = int(fields[2])
                qend = int(fields[3])
                matches = int(fields[9])
                aln_block_len = int(fields[10])
            except ValueError:
                # Skip malformed records
                continue

            identity = (matches / aln_block_len) if aln_block_len > 0 else 0.0
            left_clip = qstart
            right_clip = qlen - qend
            aln_len = aln_block_len  # use GAF column 10 (matches+mismatches+indels)

            # Keep (retain) the read UNLESS it looks plastid-like:
            # plastid-like if: both clips are small enough, identity high enough, and alignment long enough
            is_plastid_like = (
                left_clip <= clip_thr
                and right_clip <= clip_thr
                and identity > identity_thr
                and aln_len >= min_aln_len
            )
            keep = not is_plastid_like

            if keep:
                retain_out.write(qname + "\n")
            else:
                reject_out.write(qname + "\n")

            summary_writer.writerow(
                [
                    qname,
                    f"{identity:.5f}",
                    left_clip,
                    right_clip,
                    aln_len,
                    "yes" if keep else "no",
                ]
            )

    # === Close ===
    reject_out.close()
    retain_out.close()
    summary_out.close()
    print("✅ Done parsing and filtering GAF.")


if __name__ == "__main__":
    main()
