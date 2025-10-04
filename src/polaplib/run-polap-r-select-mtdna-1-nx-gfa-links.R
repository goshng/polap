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

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tidyr"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  input1 <- args[1]
  output1 <- args[2]
} else {
  s <- "Brassica_rapa"
  s <- "Vigna_radiata"
  s <- "Anthoceros_angustus"
  s <- "Spirodela_polyrhiza"
  input_dir0 <- paste0("/media/h2/goshng/figshare/", s, "/o/1/mtdna")
  input1 <- paste0(input_dir0, "/1-gfa.links.tsv")
  output1 <- paste0(input_dir0, "/1-gfa.links.edges.txt")
}


# Load the TSV file
file_path <- input1
data <- read_tsv(file_path, col_names = FALSE)

# Remove 'edge_' prefix from columns 2 and 4, then concatenate columns 2-3 and 4-5
data_processed <- data |>
  mutate(
    col1 = str_replace(X2, "edge_", "") |> paste0(X3),
    col2 = str_replace(X4, "edge_", "") |> paste0(X5)
  )

# Switch columns and signs
switch_signs <- function(x) {
  x |>
    str_replace_all("\\+", "temp") |>
    str_replace_all("-", "\\+") |>
    str_replace_all("temp", "-")
}

data_switched <- data_processed |>
  mutate(
    Source = switch_signs(col2),
    Target = switch_signs(col1)
  ) |>
  select(Source, Target)

data_processed_renamed <- data_processed |>
  mutate(Source = col1, Target = col2) |>
  select(Source, Target)

combined_data <- bind_rows(data_processed_renamed, data_switched)

# Write the processed file to a new TSV
write_tsv(combined_data, output1)
