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

################################################################################
# Copy this R template file to create a new R script that is used in polap.
# The following template file provides argument processing.
################################################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("ggplot2"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

infer_pacbio_type_from_seqkit <- function(stats_file, output_file = NULL) {
  # Detect delimiter (assumes .tsv or .csv)
  delim <- if (grepl("\\.tsv$", stats_file)) "\t" else ","

  # Read file
  df <- read_delim(stats_file, delim = delim, show_col_types = FALSE)

  # Standardize column names
  colnames(df) <- colnames(df) |>
    str_replace_all("[^A-Za-z0-9]", "_") |>
    str_replace_all("_+", "_") |>
    str_replace_all("_$", "")

  # Identify relevant columns by name
  avglen_col <- names(df)[str_detect(names(df), regex("avg.*len", ignore_case = T))][1]
  avgqual_col <- names(df)[str_detect(names(df), regex("avg.*qual", ignore_case = T))][1]

  # Convert to numeric
  df[avglen_col] <- as.numeric(df[avglen_col])
  df[avgqual_col] <- as.numeric(df[avgqual_col])

  # Classification logic: quality first, length as secondary
  df <- df %>%
    mutate(InferredType = case_when(
      !!sym(avgqual_col) >= 20 ~ "HiFi",
      !!sym(avgqual_col) < 15 ~ "CLR",
      !!sym(avgqual_col) >= 15 & !!sym(avgqual_col) < 20 &
        !!sym(avglen_col) >= 2000 ~ "Likely CLR",
      TRUE ~ "Unknown"
    )) %>%
    select(InferredType)

  # Write to file if requested
  if (!is.null(output_file)) {
    write_csv(df, output_file, col_names = FALSE)
  }

  return(df)
}

# result <- infer_pacbio_type_from_seqkit("seqkit_stats.csv")
# print(result)

# infer_pacbio_type_from_seqkit("seqkit_stats.csv", "seqkit_with_type.csv")

# args = commandArgs(trailingOnly=TRUE)
parser <- OptionParser()
parser <- add_option(parser, c("-t", "--tsv"),
  action = "store",
  help = "seqkit stats -aT TSV output filename",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = ""
)
args1 <- parse_args(parser)

if (is_null(args1$tsv)) {
  s <- "bioprojects"
  o <- "PRJNA817235-Canavalia_ensiformis"

  # input_dir0 <- file.path("/media/h2/goshng/figshare", s, o, "0")
  # test first in github/trash or github/test or somewhere at the same levels
  # of github/src.
  # Or where you have your data to work on.
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.txt")
  output1 <- file.path(input_dir0, "assembly_info_organelle_annotation_count-all.pdf")

  args1 <- parse_args(parser, args = c("--table", input1, "-o", output1))
}

infer_pacbio_type_from_seqkit(args1$tsv, args1$out)
