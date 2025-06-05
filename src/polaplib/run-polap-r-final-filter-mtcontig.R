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

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

parser <- OptionParser()
parser <- add_option(parser, c("-t", "--table"),
  action = "store",
  help = "Organelle contig annotation all table: assembly_info_organelle_annotation_count-all.txt",
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
parser <- add_option(parser, c("-l", "--length"),
  type = "integer",
  default = -1,
  help = "Maximum length of seed contigs",
  metavar = "number"
)
args1 <- parse_args(parser)

# for testing
if (is_null(args1$table)) {
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  input2 <- file.path(input_dir0, "mt.contig.name-1")
  output1 <- file.path(input_dir0, "8-mtcontig.table.tsv")

  args1 <- parse_args(parser, args = c(
    "--table", input1,
    "-m", input2,
    "-o", output1,
    "-l", 100000
  ))
}

output1 <- args1$out

# "assembly_info_organelle_annotation_count-all.txt"
x0 <- read_delim(args1$table, delim = " ", show_col_types = FALSE)
# "mt.contig.name-1"
x1 <- read_delim(args1$`mt-contig-name`, delim = " ", col_names = c("edgename"), show_col_types = FALSE)

# FIXME: we need an edge version.
xall <- x0 |>
  mutate(edgename = paste0("edge_", Edge)) |>
  relocate(edgename) |>
  arrange(Copy)

if (args1$length > 0) {
  xall <- xall |>
    filter(Length < args1$length)
}

xt <- left_join(x1, xall) |>
  filter(!is.na(MT)) |>
  arrange(desc(MT > PT), desc(MT)) |>
  select(-edgename)

# Extract the directory part
dir_path <- dirname(output1)

if (file.access(dir_path, 2) == 0) {
  cat("Directory is writable.\n")

  xt |> write_tsv(output1, col_names = FALSE)
} else {
  cat("Directory is not writable.\n")
}
