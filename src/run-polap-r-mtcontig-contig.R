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

# "$script_dir"/run-polap-mtcontig.R "${_arg_outdir}" \
# 	"${_arg_outdir}"/50-annotation/mt.contig.name \
# 	"${_arg_outdir}"/assembly_info_organelle_annotation_count.txt \
# 	--contigger \
# 	>/dev/null 2>&1
# echo "USE: assembly graph: "${_arg_outdir}"/30-contigger/graph_final.gfa"
# echo "USE: execute $ column -t "${_arg_outdir}"/assembly_info_organelle_annotation_count.txt"
# echo "INFO: edit "${_arg_outdir}"/50-annotation/mt.contig.name-1 for mtDNA contig candidates"

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))

args <- commandArgs(trailingOnly = TRUE)

# x = as_tibble(read.table(args[1]))
# copy_number_min = median(x$V6)

# consider some of these:
#   1. contig length: V2 > 1000
#   2. copy number: 1 < V6
#   3. gene count: mt > pt
# x %>% filter(V2>1000) %>% filter(V6>copy_number_min) %>% filter(mt>pt) %>% select(V1) %>% write.table(args[2],row.names=F,col.names=F,quote=F)

# consider some of these:
#   1. contig length: V2 > 1000
#   2. copy number: 1 < V6
#   3. gene count: mt > 0
# x %>% filter(V2>1000) %>% filter(V6>copy_number_min) %>% filter(mt>0) %>% select(V1) %>% write.table(args[2],row.names=F,col.names=F,quote=F)

flye_dir <- args[1]
is.contigger <- args[5]

contig.table <- paste0(flye_dir, "/contig-annotation-table.txt")
contig.depth.table <- paste0(flye_dir, "/contig-annotation-depth-table.txt")

assembly_info <- ""
if (is.contigger == "--contigger") {
  assembly_info <- paste0(flye_dir, "/30-contigger/contigs_stats.txt")
} else {
  assembly_info <- paste0(flye_dir, "/assembly_info.txt")
}

mt.gene.count <- paste0(flye_dir, "/50-annotation/mt.gene.count")
pt.gene.count <- paste0(flye_dir, "/50-annotation/pt.gene.count")

x <- as_tibble(read.table(assembly_info))
y.mt <- as_tibble(read.table(mt.gene.count)) %>% rename(mt = V2)
y.pt <- as_tibble(read.table(pt.gene.count)) %>% rename(pt = V2)
z <- x %>%
  left_join(y.mt) %>%
  left_join(y.pt) %>%
  arrange(desc(mt)) %>%
  filter(!is.na(mt)) %>%
  filter(mt > 0 | pt > 0)
z.1 <- x %>%
  left_join(y.mt) %>%
  left_join(y.pt) %>%
  arrange(V6)
z.v <- x %>%
  left_join(y.mt) %>%
  left_join(y.pt) %>%
  arrange(mt) %>%
  filter(is.na(mt))
# print("print variable: z")

mt.contig.name <- paste0(args[2], "-1")
pt.contig.name <- paste0(args[2], "-2")

copy_number_min <- median(z.1$V6)

# z %>% filter(V2>1000) %>% filter(V6>copy_number_min) %>% filter(mt>pt) %>% select(V1) %>% write.table(mt.contig.name,row.names=F,col.names=F,quote=F)
# z %>% filter(V2>1000) %>% filter(V6>copy_number_min) %>% filter(mt<pt) %>% select(V1) %>% write.table(pt.contig.name,row.names=F,col.names=F,quote=F)

z %>%
  arrange(mt <= pt) %>%
  relocate(V9, .after = last_col()) %>%
  rename(Contig = V1, Length = V2, Copy = V6, MT = mt, PT = pt, Edge = V9) %>%
  write.table(args[3], row.names = F, quote = F)
z %>%
  arrange(mt <= pt) %>%
  relocate(V9, .after = last_col()) %>%
  filter(mt > pt) %>%
  select(V1, V2, V6, mt, pt, V9) %>%
  rename(Contig = V1, Length = V2, Copy = V6, MT = mt, PT = pt, Edge = V9) %>%
  write.table(contig.table, row.names = F, quote = F)
z %>%
  arrange(mt <= pt) %>%
  relocate(V9, .after = last_col()) %>%
  filter(mt > pt) %>%
  select(V1, V2, V3, V6, mt, pt, V9) %>%
  rename(Contig = V1, Length = V2, Depth = V3, Copy = V6, MT = mt, PT = pt, Edge = V9) %>%
  write.table(contig.depth.table, row.names = F, quote = F)
z.1 %>%
  rename(Contig = V1, Length = V2, Copy = V6, MT = mt, PT = pt, Edge = V9) %>%
  relocate(Edge, .after = last_col()) %>%
  write.table(args[4], row.names = F, quote = F)
