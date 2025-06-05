#!/usr/bin/env Rscript

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

# name: converts BLAST output file to a BED format file.
#
# synopsis:
# run-polap-genes.R $MTAABLAST $MTAABLAST.bed >/dev/null 2>&1
# input: MTAABLAST
# output: BED
#
# requirement:
# Create a tblastn output file: MTAABLAST
# tblastn -query $MTAA -db $CONTIGDB -out $MTAABLAST -evalue 1e-30 \
# 	-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle salltitles"

suppressPackageStartupMessages(library("dplyr"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

args <- commandArgs(trailingOnly = TRUE)

x <- as_tibble(read.table(args[1]))
y1 <- x %>%
  select(V1, V2, V9, V10) %>%
  filter(V9 < V10) %>%
  rename(X1 = V2, X2 = V9, X3 = V10, X4 = V1) %>%
  relocate(X4, .after = last_col())
y2 <- x %>%
  select(V1, V2, V9, V10) %>%
  filter(V9 > V10) %>%
  rename(X1 = V2, X2 = V10, X3 = V9, X4 = V1) %>%
  relocate(X4, .after = last_col())
z <- bind_rows(y1, y2)

write.table(z, args[2], row.names = F, col.names = F, quote = F, sep = "\t")
