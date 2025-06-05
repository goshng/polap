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
parser <- add_option(parser, c("--unit-annotation-edge"),
  action = "store_true",
  default = FALSE, help = "true for edge, false for contig"
)
parser <- add_option(parser, c("-m", "--mitochondrial"),
  action = "store_true",
  default = TRUE, help = "Mitochondrial genome assembly"
)
parser <- add_option(parser, c("-p", "--plastid"),
  action = "store_false",
  dest = "mitochondrial", help = "Plastid genome assembly"
)
parser <- add_option(parser, c("-g", "--gfa"),
  action = "store",
  help = "GFA sequence part",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("--edge"),
  action = "store",
  help = "Depth range",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output contig seeds filename"
)
args1 <- parse_args(parser)

if (is_null(args1$gfa)) {
  s <- "bioprojects"
  o <- "PRJNA817235-Canavalia_ensiformis"

  # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "4-gfa.links.tsv")
  input2 <- file.path(input_dir0, "1-preselection.by.gene.density.txt")
  # input2 <- file.path(input_dir0, "1-preselection.by.depth.mixture.txt")
  output.base <- file.path(input_dir0, "4-gfa.links")
  args1 <- parse_args(parser, args = c("--gfa", input1, "--edge", input2, "-o", output.base))
}

output1 <- paste0(args1$out, ".number.txt")
output2 <- paste0(args1$out, ".order.txt")
output3 <- paste0(args1$out, ".contig.txt")
output4 <- paste0(args1$out, ".contig.na.txt")

# x1 <- read_tsv(args1$edge, col_names = c("string", "depth"))
#

if (args1$`unit-annotation-edge` == TRUE) {
  x1 <- read_tsv(args1$edge, col_names = c("string", "length", "depth", "V4", "V5", "copy", "V7", "V8", "mt", "pt", "edge", "density"), show_col_types = FALSE)
} else {
  x1 <- read_tsv(args1$edge, col_names = c("string", "contig", "length", "depth", "copy", "mt", "pt", "edge", "density"), show_col_types = FALSE)
}

# Step 1: Read the TSV file (assuming it has no header and two columns)
# Adjust the file path as necessary
df <- read_tsv(args1$gfa, col_names = c("string1", "string2"), show_col_types = FALSE)

# Step 2: Combine both columns to identify all unique strings
unique_strings <- df |>
  dplyr::select(string1, string2) |>
  unlist() |>
  unique() |>
  sort()

# Step 3: Create a mapping table: each string corresponds to a unique non-negative integer
mapping_table <- tibble(
  int_value = seq(0, length(unique_strings) - 1),
  string = unique_strings
)

vector_mapped <- x1 |>
  left_join(mapping_table, by = "string") |>
  filter(!is.na(int_value)) |>
  dplyr::select(int_value)

vector_mapped_na <- x1 |>
  left_join(mapping_table, by = "string") |>
  filter(is.na(int_value)) |>
  dplyr::select(string)

# Step 4: Use left_join to map each string to its corresponding integer
df_mapped <- df |>
  left_join(mapping_table, by = c("string1" = "string")) |>
  rename(int1 = int_value) |>
  left_join(mapping_table, by = c("string2" = "string")) |>
  rename(int2 = int_value) |>
  dplyr::select(int1, int2)



# Step 5: Now `df_mapped` contains integer values corresponding to each string
# and `mapping_table` shows the string-to-integer mapping

# Show the first few rows of the result
df_mapped |>
  write_tsv(output1, col_names = FALSE)

# Show the mapping table
mapping_table |>
  write_tsv(output2, col_names = FALSE)

vector_mapped |>
  write_tsv(output3, col_names = FALSE)

vector_mapped_na |>
  write_tsv(output4, col_names = FALSE)
