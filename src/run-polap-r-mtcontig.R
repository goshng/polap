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

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))

parser <- OptionParser()

parser <- add_option(parser, c("--flyeout-edges-stats"),
                     action = "store",
                     help = "Input1 30-contigger/edges_stats.txt"
)
parser <- add_option(parser, c("--mt-gene-count"),
                     action = "store",
                     help = "Input2 50-annotation/mt.gene.count"
)
parser <- add_option(parser, c("--pt-gene-count"),
                     action = "store",
                     help = "Input3 50-annotation/pt.gene.count"
)
parser <- add_option(parser, c("--out-annotation"),
                     action = "store",
                     help = "Output1 assembly_info_organelle_annotation_count.txt"
)
parser <- add_option(parser, c("--out-annotation-all"),
                     action = "store",
                     help = "Output2 assembly_info_organelle_annotation_count-all.txt"
)
parser <- add_option(parser, c("--out-annotation-table"),
                     action = "store",
                     help = "Output3 contig-annotation-table.txt"
)
parser <- add_option(parser, c("--out-annotation-depth-table"),
                     action = "store",
                     help = "Output4 contig-annotation-depth-table.txt"
)
parser <- add_option(parser, c("-c", "--contigger"),
                     action = "store_true",
                     default = FALSE,
                     help = "Use 30-contigger flye output"
)
args1 <- parse_args(parser)

if (is_null(args1$`flyeout-edges-stats`)) {
  s <- "bioprojects"
  o <- "PRJNA817235-Canavalia_ensiformis"

  # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  # test first in github/trash or github/test or somewhere at the same levels
  # of github/src.
  # Or where you have your data to work on.
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "edges_stats.txt")
  input2 <- file.path(input_dir0, "mt.gene.count")
  input3 <- file.path(input_dir0, "pt.gene.count")
  output1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count.txt")
  output2 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  output3 <- file.path(input_dir0, "contig-annotation-table.txt")
  output4 <- file.path(input_dir0, "contig-annotation-depth-table.txt")
  
  args1 <- parse_args(parser, args = c("--flyeout-edges-stats", input1, 
                                       "--mt-gene-count", input2, 
                                       "--pt-gene-count", input3, 
                                       "--out-annotation", output1, 
                                       "--out-annotation-all", output2, 
                                       "--out-annotation-table", output3, 
                                       "--out-annotation-depth-table", output4, 
                                       "--contigger"))
}

# flye_dir <- args[1]
# is.contigger <- args[5]

contig.table <- args1$`out-annotation-table`
contig.depth.table <- args1$`out-annotation-depth-table`

assembly_info <- args1$`flyeout-edges-stats`
# if (args1$contigger == TRUE) {
#   assembly_info <- paste0(flye_dir, "/30-contigger/contigs_stats.txt")
# } else {
#   assembly_info <- paste0(flye_dir, "/assembly_info.txt")
# }

mt.gene.count <- args1$`mt-gene-count`
pt.gene.count <- args1$`pt-gene-count`

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

# mt.contig.name <- paste0(args[2], "-1")
# pt.contig.name <- paste0(args[2], "-2")

copy_number_min <- median(z.1$V6)

# z %>% filter(V2>1000) %>% filter(V6>copy_number_min) %>% filter(mt>pt) %>% select(V1) %>% write.table(mt.contig.name,row.names=F,col.names=F,quote=F)
# z %>% filter(V2>1000) %>% filter(V6>copy_number_min) %>% filter(mt<pt) %>% select(V1) %>% write.table(pt.contig.name,row.names=F,col.names=F,quote=F)

z %>%
  arrange(mt <= pt) %>%
  relocate(V9, .after = last_col()) %>%
  rename(Contig = V1, Length = V2, Copy = V6, MT = mt, PT = pt, Edge = V9) %>%
  write.table(args1$`out-annotation`, row.names = F, quote = F)

z %>%
  arrange(mt <= pt) %>%
  relocate(V9, .after = last_col()) %>%
  filter(mt > pt) %>%
  select(V1, V2, V6, mt, pt, V9) %>%
  rename(Contig = V1, Length = V2, Copy = V6, MT = mt, PT = pt, Edge = V9) %>%
  write.table(args1$`out-annotation-table`, row.names = F, quote = F)

z %>%
  arrange(mt <= pt) %>%
  relocate(V9, .after = last_col()) %>%
  filter(mt > pt) %>%
  select(V1, V2, V3, V6, mt, pt, V9) %>%
  rename(Contig = V1, Length = V2, Depth = V3, Copy = V6, MT = mt, PT = pt, Edge = V9) %>%
  write.table(args1$`out-annotation-depth-table`, row.names = F, quote = F)

z.1 %>%
  rename(Contig = V1, Length = V2, Copy = V6, MT = mt, PT = pt, Edge = V9) %>%
  relocate(Edge, .after = last_col()) %>%
  write.table(args1$`out-annotation-all`, row.names = F, quote = F)
