#!/usr/bin/env python3
import sys
import csv

# === Inputs ===
gaf_path = sys.argv[1]  # input GAF file
reject_path = sys.argv[2]  # output: rejected read names
retain_path = sys.argv[3]  # output: retained read names
summary_path = sys.argv[4]  # output: summary.tsv
identity_thresh = float(sys.argv[5])  # e.g., 0.95
clip_thresh = int(sys.argv[6])  # e.g., 100

# === Open files ===
reject_out = open(reject_path, "w")
retain_out = open(retain_path, "w")
summary_out = open(summary_path, "w", newline="")
summary_writer = csv.writer(summary_out, delimiter="\t")
summary_writer.writerow(["read_id", "identity", "left_clip", "right_clip", "keep"])

with open(gaf_path) as f:
    for line in f:
        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")
        if len(fields) < 12:
            continue  # not a valid alignment

        qname = fields[0]
        qlen = int(fields[1])
        qstart = int(fields[2])
        qend = int(fields[3])
        matches = int(fields[9])
        aln_block_len = int(fields[10])

        identity = matches / aln_block_len if aln_block_len > 0 else 0

        left_clip = qstart
        right_clip = qlen - qend

        keep = not (
            left_clip <= clip_thresh
            and right_clip <= clip_thresh
            and identity > identity_thresh
        )

        if keep:
            retain_out.write(qname + "\n")
        else:
            reject_out.write(qname + "\n")

        summary_writer.writerow(
            [qname, f"{identity:.5f}", left_clip, right_clip, "yes" if keep else "no"]
        )

# === Close ===
reject_out.close()
retain_out.close()
summary_out.close()
print("✅ Done parsing and filtering GAF.")
