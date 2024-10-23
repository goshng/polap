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
parser <- add_option(parser, c("-t", "--table"),
                     action = "store",
                     help = "Organelle annotation all table",
                     metavar = "<FILE>"
)
parser <- add_option(parser, c("-a", "--annotation"),
                     action = "store",
                     help = "Organelle annotation depth table",
                     metavar = "<FILE>"
)
parser <- add_option(parser, c("-m", "--mtcontigname"),
                     action = "store",
                     help = "Organelle contigs names",
                     metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output contig seeds filename"
)
args1 <- parse_args(parser)

if (is_null(args1$table)) {
  s <- "bioprojects"
  o <- "PRJNA817235-Canavalia_ensiformis"

  # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  # test first in github/trash or github/test or somewhere at the same levels
  # of github/src.
  # Or where you have your data to work on.
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  input2 <- file.path(input_dir0, "contig-annotation-depth-table.txt")
  input3 <- file.path(input_dir0, "mt.contig.name-9")
  output1 <- file.path(input_dir0, "contig-annotation-depth-table-seeds.txt")

  args1 <- parse_args(parser, args = c("--table", input1, "-a", input2, "-m", input3, "-o", output1))
}

# output1 <- paste0(args1$out, ".table.tsv")
output1 <- paste0(args1$out)

table2 <- read_delim(args1$table, delim = " ", show_col_types = FALSE)
table1 <- read_delim(args1$annotation, delim = " ", show_col_types = FALSE)

# Load table1 and table2
#table1 <- read_tsv("path/to/contig-annotation-depth-table.txt")
#table2 <- read_tsv("path/to/assembly_info_organelle_annotation_count-all.txt", col_names = FALSE)

# Load the edge names set
edges <- read_lines(args1$mtcontigname)

# Process table2 by removing columns V4, V5, V7, V8 and renaming V3 as Depth
table2_processed <- table2 %>%
  select(Contig, Length, V3, Copy, MT, PT, Edge) %>%
  rename(Depth = V3)

# Check if each Contig in table1 exists in the set of edge names
table1 <- table1 %>%
  mutate(Seed = ifelse(Contig %in% edges, "A", "X"))

# Identify edges not present in table1
missing_edges <- setdiff(edges, table1$Contig)

# For the missing edges, create new rows from table2 and set Seed to "G"
new_rows <- table2_processed %>%
  filter(Contig %in% missing_edges) %>%
  mutate(Seed = "G")

# Append the new rows to table1
final_table <- bind_rows(table1, new_rows) %>%
  select(-Edge)

# Save the final table as a TSV file
write_tsv(final_table, output1)

