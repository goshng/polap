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

# name: select long reads
#
# synopsis:
# 	run-polap-pairs.R mt.contig.name-1 contig.tab ${MTSEEDSDIR} $SINGLE_MIN $SINGLE_MIN
#
# requirement: executes Flye
# flye --nano-raw "$LR3K" --out-dir "${_arg_outdir}" \
# 	--threads "${_arg_threads}" \
# 	--stop-after contigger \
# 	--asm-coverage 30 \
# 	--genome-size "$EXPECTED_GENOME_SIZE"
#
# input:
#   1. mt.contig.name-1 - contig or edge numbers
#   2. contig.tab - minimap2 output modified
#   3. seeds directory for output
#   4. pair minimum length V11: 3000 for MT, 1000 for PT
#   5. bridge minimum length V7: depends: 3000 or 5000
#   6. single minimum length V11: 3000 for MT, 0 or 1000 for PT
#
# output:
#   single.names and <contig1-contig2.name> files in the output directory
#
# MTSEEDSDIR="${_arg_outdir}"/60-mt-${STEP4}/o${MR}/seeds
# 	run-polap-pairs.R mt.contig.name-1 contig.tab ${MTSEEDSDIR} $PAIR_MIN $BRIDGE_MIN $SINGLE_MIN
# 	run-polap-pairs.R mt.contig.name-1 contig.tab ${MTSEEDSDIR} $PAIR_MIN $BRIDGE_MIN

suppressWarnings(suppressPackageStartupMessages(library("optparse")))
suppressWarnings(suppressPackageStartupMessages(library("dplyr")))
suppressWarnings(suppressPackageStartupMessages(library("readr")))
suppressWarnings(suppressPackageStartupMessages(library("purrr")))
suppressWarnings(suppressPackageStartupMessages(library("tidyr")))
suppressWarnings(suppressPackageStartupMessages(library("ggplot2")))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

parser <- OptionParser()

clean_string <- function(input_string) {
  # Replace newlines with a space
  cleaned_string <- gsub("\\n", " ", input_string)

  # Replace multiple spaces with a single space
  cleaned_string <- gsub("\\s+", " ", cleaned_string)

  # Trim leading and trailing spaces
  cleaned_string <- trimws(cleaned_string)

  return(cleaned_string)
}

parser <- add_option(parser, c("-t", "--table"),
  type = "character",
  default = "contig1.tab",
  action = "store",
  metavar = "<FILE>",
  help = "Input minimap2 PAF-like tabular format file [default %default]"
)
parser <- add_option(parser, c("-o", "--out"),
  type = "character",
  default = "ptgaul.contig1.names.txt",
  action = "store",
  metavar = "<FILE>",
  help = "Output read name file [default %default]"
)
parser <- add_option(parser, c("-w", "--single-min"),
  type = "integer",
  default = 3000,
  metavar = "<integer>",
  help = "Minimum length of single-mapping alignment [default %default]"
)

parser <- add_option(parser, c("--ratio"),
  type = "double",
  default = 0.7,
  help = "Minimum of V10(match)/V11(base) [default %default]",
  metavar = "<real>"
)

long_help_text <- "Use only the forward-strand reads
[default %default]"
parser <- add_option(parser, c("--forward"),
  type = "logical",
  default = FALSE,
  action = "store_true",
  metavar = "BOOLEAN",
  help = paste(strwrap(clean_string(long_help_text), width = 80),
    collapse = "\n\t\t"
  )
)

long_help_text <- "Use only the reverse-strand reads
[default %default]"
parser <- add_option(parser, c("--reverse"),
  type = "logical",
  default = FALSE,
  action = "store_true",
  metavar = "BOOLEAN",
  help = paste(strwrap(clean_string(long_help_text), width = 80),
    collapse = "\n\t\t"
  )
)

args1 <- parse_args(parser)

if (is_null(args1$table)) {
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "contig1.tab")
  output1 <- file.path(input_dir0, "ptgaul.contig1.names.txt")

  args1 <- parse_args(parser, args = c(
    "--table", input1,
    "-o", output1,
    "-w", 10000,
    "--forward"
  ))
}

# Check for exclusive options logic (only one can be TRUE)
if (args1$forward && args1$reverse) {
  stop("only one of --forward or --reverse can be selected.")
}

# https://lh3.github.io/minimap2/minimap2.html
# Li 2008 - "OUTPUT FORMAT: Minimap2 outputs mapping positions in the
# Pairwise mApping Format (PAF) by default. PAF is a TAB-delimited text format
# with each line consisting of at least 12 fields as are described in the
# following table:
#
# Col	Type	Description
# 1	string	Query sequence name
# 2	int	Query sequence length
# 3	int	Query start coordinate (0-based)
# 4	int	Query end coordinate (0-based)
# 5	char	‘+’ if query/target on the same strand; ‘-’ if opposite
# 6	string	Target sequence name
# 7	int	Target sequence length
# 8	int	Target start coordinate on the original strand
# 9	int	Target end coordinate on the original strand
# 10	int	Number of matching bases in the mapping
# 11	int	Number bases, including gaps, in the mapping
# 12	int	Mapping quality (0-255 with 255 for missing)
#
# Read in the data with assigned column names
data <- read_tsv(args1$table,
  show_col_types = FALSE,
  col_names = c(
    "rname",
    "rlen",
    "rstart",
    "rend",
    "strand",
    "cname",
    "clen",
    "cstart",
    "cend",
    "match",
    "base"
  )
)

if (args1$forward) {
  data <- data |> filter(strand == "+")
} else if (args1$reverse) {
  data <- data |> filter(strand == "-")
}

single_min <- args1$`single-min`

# 0. Check preconditions
# rstart is not greater than rend.
# match is not greater than base.
result <- data |>
  summarize(all_rows_meet_condition = all(rstart <= rend & match <= base))
stopifnot(result$all_rows_meet_condition)
stopifnot(single_min >= 0)

# 1. ptGAUL: we use ptGAUL for the case of single edge reference
ptgaul_file <- args1$out
ptgaul_mapped_reads <- data |>
  filter(match / base > args1$ratio, base > single_min) |>
  select(rname) |>
  distinct(rname)

ptgaul_mapped_reads |>
  write.table(
    ptgaul_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

# ptgaul_option_file <- file.path(mtdir, ptgaul_option_base)
# print(ptgaul_file)
# print(ptgaul_option_file)
# file.copy(ptgaul_file, ptgaul_option_file)
# file.copy(ptgaul_file, ptgaul_option_file, showWarnings = FALSE)
