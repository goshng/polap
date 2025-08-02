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

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
})

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

# Define options
option_list <- list(
  make_option(c("-p", "--min-identity"),
    type = "double", default = 0.80,
    help = "Minimum identity threshold [default %default]"
  ),
  make_option(c("-w", "--min-length"),
    type = "integer", default = 1000,
    help = "Minimum alignment length [default %default]"
  )
)

# Parse options and arguments
parser <- OptionParser(
  usage = "%prog [options] <input.paf> <output.txt>",
  option_list = option_list
)
args <- parse_args(parser, positional_arguments = 2)

# if (is_null(args1$args[1])) {
#   s <- "bioprojects"
#   o <- "PRJNA817235-Canavalia_ensiformis"

#   # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
#   # test first in github/trash or github/test or somewhere at the same levels
#   # of github/src.
#   # Or where you have your data to work on.
#   input_dir0 <- file.path(".")
#   input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
#   output1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.pdf")

#   args1 <- parse_args(parser, args = c("--table", input1, "-o", output1))
# }

input_paf <- args$args[1]
output_txt <- args$args[2]
min_identity <- args$options$`min-identity`
min_length <- args$options$`min-length`

# Read PAF file
paf <- read_tsv(input_paf, col_names = FALSE, comment = "@", show_col_types = FALSE)

# Assign column names for alignment fields
colnames(paf)[1:11] <- c(
  "qname", "qlen", "qstart", "qend",
  "strand", "tname", "tlen", "tstart", "tend",
  "nmatch", "alen"
)

# Filter reads by alignment length and identity
paf_filtered <- paf %>%
  mutate(identity = nmatch / alen) %>%
  filter(alen >= min_length, identity >= min_identity)

# Write unique read names
write_lines(unique(paf_filtered$qname), output_txt)
