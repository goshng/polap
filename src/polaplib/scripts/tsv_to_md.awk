#!/usr/bin/env awk
# scripts/tsv_to_md.awk
# Version: v0.1.0
# SPDX-License-Identifier: GPL-3.0-or-later
BEGIN { FS="\t"; OFS="|" }
NR==1 {
  printf("|")
  for (i=1;i<=NF;i++) printf(" %s |",$i)
  printf("\n|")
  for (i=1;i<=NF;i++) printf(" --- |")
  printf("\n")
  next
}
{
  printf("|")
  for (i=1;i<=NF;i++) printf(" %s |",$i)
  printf("\n")
}
