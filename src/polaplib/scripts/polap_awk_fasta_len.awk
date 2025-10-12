#!/usr/bin/env awk -f
# Version: v0.8.0
# Sum total length of sequences in a FASTA
# Usage: awk -f polap_awk_fasta_len.awk mt.fasta
/^>/ {next}
{ gsub(/\r/,""); L+=length($0) }
END { print L+0 }

