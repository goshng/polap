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

# polaplib/polap-r-final-seed-mtcontig.R

# TODO: document it later
#
# Subcommand: seeds
#
# Used by:
# function _polap_seeds_final-mtcontig {
# _polap_disassemble_seeds_final-seeds-mtcontig() {
#
# See Also:
# polap-r-final-seed-mtcontig.R
# function _polap_disassemble_seeds_final-seeds-mtcontig : polap-function-disassemble-seeds.sh
# function _polap_seeds_final-seeds-mtcontig : run-polap-function-seeds.sh
# subcommand seeds uses this script indirectly.
#

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("tidyr"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

parser <- OptionParser()
parser <- add_option(parser, c("-t", "--table"),
  action = "store",
  help = "Organelle contig annotation all table: assembly_info_organelle_annotation_count-all.txt",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-a", "--annotation"),
  action = "store",
  help = "Organelle contig annotation depth table: contig-annotation-depth-table.txt",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-m", "--mt-contig-name"),
  action = "store",
  help = "Organelle seed contigs: mt.contig.name-1",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output contig table: mtcontig.table.tsv",
  metavar = "<FILE>"
)

args1 <- parse_args(parser)

# for testing
if (is_null(args1$table)) {
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "jvalidus/disassemble/8/assembly_info_organelle_annotation_count-all.txt")
  input2 <- file.path(input_dir0, "jvalidus/disassemble/8/51-mtcontigs/2/7-mt.contig.name.txt")
  input3 <- file.path(input_dir0, "jvalidus/disassemble/8/contig-annotation-depth-table.txt")
  output1 <- file.path(input_dir0, "8-mtcontig.table.tsv")
  output2 <- file.path(input_dir0, "jvalidus/disassemble/8/51-mtcontigs/2/8-mtcontig-annotation-table-seed.txt")

  args1 <- parse_args(parser, args = c(
    "--table", input1,
    "-m", input2,
    "-a", input3,
    "-o", output2
  ))
}

output2 <- paste0(args1$out)

# table2: "assembly_info_organelle_annotation_count-all.txt"
table2 <- read_delim(args1$table, delim = " ", show_col_types = FALSE)
# table1: contig-annotation-depth-table.txt
table1 <- read_delim(args1$annotation, delim = " ", show_col_types = FALSE)
# edges: mt.conting.name-1
edges <- read_lines(args1$`mt-contig-name`)

# Process table2 by removing columns V4, V5, V7, V8 and renaming Depth as Depth
table2_processed <- table2 %>%
  select(Contig, Length, Depth, Copy, MT, PT, Edge)

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
write_tsv(final_table, output2)
