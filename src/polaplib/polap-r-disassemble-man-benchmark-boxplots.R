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
# This creates the three different kinds of  figures of processing time and 
# memory performance
# of different ptDNA assembly pipelines: ptGAUL, GetOrganelle, PMAT, Oatk, and
# TIPPo.
# It reads a table with the results of the benchmarking and creates boxplots for
# the time and memory usage of each method. We have three different types of 
# plots:
# 1. Time for the assembly of the ptDNA
# 2. Memory usage during the assembly of the ptDNA
# 3. Time for the assembly of the ptDNA with NextDenovo
#
# The input table is expected to have the following columns:
# - `_species`: the species name
# - `_total_hours_*`: the total time in hours for each method
# - `_memory_gb_*`: the memory usage in GB for each method
# The output is a PDF file with the boxplots for each method.
#
# The methods are ordered and labeled as follows:
# - "genomesize" = "G"
# - "msbwt" = "MSBWT"
# - "msbwt_polish" = "FMLRC"
# - "nextdenovo_polish" = "NextDenovo"
# - "getorganelle" = "GetOrganelle"
# - "ptgaul" = "ptGAUL"
# - "polap_disassemble_default" = "Polap"
# - "oatk_nextdenovo_30" = "Oatk-30"
# - "oatk_nextdenovo_20" = "Oatk-20"
# - "tippo_nextdenovo_onthq" = "TIPPo-hq"
# - "tippo_nextdenovo_ont" = "TIPPo-ont"
# - "pmat_nextdenovo_0.1" = "PMAT-0.1"
# - "pmat_nextdenovo_1.0" = "PMAT-1.0"
#
# Note: we convert the time from hours to minutes for the plots: e.g.,
# 3h or 15m will be converted to 180 minutes or 15 minutes, respectively.
# We do not convert 1h 30m to 90 minutes. So, the time input format must be
# either "3h" or "15m".
#
# Example:
# polap-data-cflye man figure-benchmark some 2 time
# polap-data-cflye man figure-benchmark some 2 memory
# polap-data-cflye man figure-benchmark some 2 time-nextdenovo
#
# Rscript polap-r-disassemble-man-benchmark-boxplots.R \
#     --table input/table-some-2.tsv \
#     --type time \
#     -o figure-time-some-2.pdf
#
# Rscript polap-r-disassemble-man-benchmark-boxplots.R \
#     --table input/table-some-2.tsv \
#     --type memory \
#     -o figure-memory-some-2.pdf
#
# Rscript polap-r-disassemble-man-benchmark-boxplots.R \
#     --table input/table-some-2.tsv \
#     --type time-nextdenovo \
#     -o figure-time-nextdenovo-some-2.pdf
#
# Tips:
# We use pivot_longer to transform the wide format of the input table.
# We use map_dbl to convert the time format to apply functions to each row.
#
# Check: 2025-06-20
################################################################################

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

# We use optparse to parse the command line arguments.
# The arguments are:
# -t or --table: the input table with the benchmarking results (required)
# -o or --out: the output PDF filename (required)
# -d or --type: the type of benchmarking to plot (time, memory, or time-nextdenovo)
# If no arguments are provided, we use a default input file for testing
# purposes.
parser <- OptionParser()
parser <- add_option(parser, c("-t", "--table"),
  action = "store",
  type = "character",
  default = NULL,
  help = "Benchmarking table [required]",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  type = "character",
  default = NULL,
  help = "Output PDF filename [required]",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-d", "--type"),
  action = "store",
  type = "character",
  default = "time",
  help = "Benchmarking type: time, memory, or time-nextdenovo",
  metavar = "<TYPE>"
)
args1 <- parse_args(parser)

# If no table is provided, use a default input file for testing
# This is useful for development and testing purposes.
# If you run this script without any arguments, 
# it will use the default input file.
if (is_null(args1$table)) {
  
  input_dir0 <- file.path("input")
  input_dir1 <- file.path("output")
  input1 <- file.path(input_dir0, "table-some-2.tsv")
  
  # for time test
  input2 <- "time"
  output1 <- file.path(input_dir1, "figure-time-some-2.pdf")

  # for time-nextdenovo test
  input2 <- "time-nextdenovo"
  output1 <- file.path(input_dir1, "figure-time-nextdenovo-some-2.pdf")

  # for memory test
  input2 <- "memory"
  output1 <- file.path(input_dir1, "figure-memory-some-2.pdf")

  # for time test
  input2 <- "time"
  output1 <- file.path(input_dir1, "figure-time-some-2.pdf")
  
  args1 <- parse_args(parser, args = c(
    "--table", input1,
    "-o", output1,
    "--type", input2
  ))
}

# Check required options
required_opts <- c("table", "out")
missing <- required_opts[sapply(required_opts, function(x) is.null(args1[[x]]))]

if (length(missing) > 0) {
  print_help(parser)
  stop(paste("Missing required option(s):", paste(missing, collapse=", ")), call.=FALSE)
}

# Check if the input file exists
if (!file.exists(args1$table)) {
  stop(paste("Input file does not exist:", args1$table))
}

# Check if the output file is provided
if (is.null(args1$out)) {
  stop("Output file must be specified with -o or --out")
}

################################################################################
# Main script logic
#
output1 <- paste0(args1$out)

# Load data
df <- read_tsv(args1$table, show_col_types = FALSE)

# We have two types of plots: time and memory.
# If the type is "time-nextdenovo", we will use the same logic as for "time",
# but we will filter the methods differently.
# We have three types of plots using this if-statement.

################################################################################
# If the type is "memory", we will plot the memory usage.
#
if (args1$type == "memory") {
  # Define method order and labels
  method_order <- c(
    "genomesize",
    "msbwt",
    "msbwt_polish",
    "nextdenovo_polish",
    "getorganelle",
    "ptgaul",
    "polap_disassemble_default",
    "polap_disassemble_simple",
    "oatk_nextdenovo_30",
    "oatk_nextdenovo_20",
    "tippo_nextdenovo_onthq",
    "tippo_nextdenovo_ont",
    "pmat_nextdenovo_0.1",
    "pmat_nextdenovo_1.0"
  )

  method_labels <- c(
    "genomesize" = "G",
    "msbwt" = "MSBWT",
    "msbwt_polish" = "FMLRC",
    "nextdenovo_polish" = "NextDenovo",
    "getorganelle" = "GetOrganelle",
    "ptgaul" = "ptGAUL",
    "polap_disassemble_default" = "Polap",
    "polap_disassemble_simple" = "Polap-ptGAUL",
    "oatk_nextdenovo_30" = "Oatk-30",
    "oatk_nextdenovo_20" = "Oatk-20",
    "tippo_nextdenovo_onthq" = "TIPPo-hq",
    "tippo_nextdenovo_ont" = "TIPPo-ont",
    "pmat_nextdenovo_0.1" = "PMAT-0.1",
    "pmat_nextdenovo_1.0" = "PMAT-1.0"
  )

  # We will remove the following methods from the memory plot
  methods_to_remove <- c(
    "polap_disassemble_simple",
    "polap_disassemble_check_default",
    "polap_disassemble_check_simple"
  )

  # Process memory data
  df_memory <- df %>%
    select(`_species`, starts_with("_memory_gb_")) %>%
    pivot_longer(
      cols = starts_with("_memory_gb_"),
      names_to = "method",
      values_to = "memory_gb"
    ) %>%
    mutate(method = str_remove(method, "_memory_gb_")) %>%
    filter(!method %in% methods_to_remove) %>%
    drop_na(memory_gb) %>%
    filter(method %in% method_order) %>%
    mutate(
      method_label = factor(method_labels[method], levels = method_labels[method_order])
    )

  # Plot
  p1 <- ggplot(df_memory, aes(x = method_label, y = memory_gb)) +
    geom_boxplot(fill = "gray80", color = "black") + # Light gray fill
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(
      # title = "Memory Usage by Method",
      x = "Method",
      y = "Memory (GB)"
    )

  # Save the plot as a PDF
  ggsave(output1, plot = p1, device = "pdf", width = 8, height = 6)
} else {
  ################################################################################
  # If the type is "time", we will plot the total time in minutes for each method.
  #
  
  # Helper to convert time to minutes
  convert_time_to_minutes <- function(x) {
    if (is.na(x)) {
      return(NA_real_)
    }
    x <- str_trim(x)
    if (str_detect(x, "h$")) {
      as.numeric(str_remove(x, "h$")) * 60
    } else if (str_detect(x, "m$")) {
      as.numeric(str_remove(x, "m$"))
    } else {
      NA_real_
    }
  }

  # Define method order and labels
  method_order <- c(
    "genomesize",
		"msbwt",
		"msbwt_polish",
		"nextdenovo_polish",
    "getorganelle",
		"ptgaul",
		"polap_disassemble_default",
		"polap_disassemble_simple",
    "oatk_nextdenovo_30",
		"oatk_nextdenovo_20",
    "tippo_nextdenovo_onthq",
		"tippo_nextdenovo_ont",
    "pmat_nextdenovo_0.1",
		"pmat_nextdenovo_1.0"
  )

  method_labels <- c(
    "genomesize" = "G",
    "msbwt" = "MSBWT",
    "msbwt_polish" = "FMLRC",
    "nextdenovo_polish" = "NextDenovo",
    "getorganelle" = "GetOrganelle",
    "ptgaul" = "ptGAUL",
    "polap_disassemble_default" = "Polap",
    "polap_disassemble_simple" = "Polap-ptGAUL",
    "oatk_nextdenovo_30" = "Oatk-30",
    "oatk_nextdenovo_20" = "Oatk-20",
    "tippo_nextdenovo_onthq" = "TIPPo-hq",
    "tippo_nextdenovo_ont" = "TIPPo-ont",
    "pmat_nextdenovo_0.1" = "PMAT-0.1",
    "pmat_nextdenovo_1.0" = "PMAT-1.0"
  )

  if (args1$type == "time") {
    methods_to_remove <- c(
      "nextdenovo_polish",
      "polap_disassemble_simple",
      "polap_disassemble_check_default",
      "polap_disassemble_check_simple"
    )
  } else {
    methods_to_remove <- c(
      "polap_disassemble_simple",
      "polap_disassemble_check_default",
      "polap_disassemble_check_simple"
    )
  }

  # Process and filter time data
  df_time <- df %>%
    select(`_species`, starts_with("_total_hours_")) %>%
    pivot_longer(
      cols = starts_with("_total_hours_"),
      names_to = "method",
      values_to = "time_raw"
    ) %>%
    mutate(
      method = str_remove(method, "_total_hours_"),
      total_minutes = map_dbl(time_raw, convert_time_to_minutes)
    ) %>%
    filter(!method %in% methods_to_remove) %>%
    drop_na(total_minutes) %>%
    filter(method %in% method_order) %>%
    mutate(
      method_label = factor(method_labels[method], levels = method_labels[method_order])
    )

  # Plot
  p1 <- ggplot(df_time, aes(x = method_label, y = total_minutes)) +
    geom_boxplot(fill = "gray80", color = "black") + # Light gray fill
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(
      # title = "Total Time by Method",
      x = "Method",
      y = "Time (minutes)"
    )

  # Save the plot as a PDF
  ggsave(output1, plot = p1, device = "pdf", width = 8, height = 6)
}
