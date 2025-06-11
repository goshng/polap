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
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("ggplot2"))

debug <- Sys.getenv("_POLAP_DEBUG", unset = "0")

# args = commandArgs(trailingOnly=TRUE)
parser <- OptionParser()
parser <- add_option(parser, c("-t", "--table"),
  action = "store",
  help = "Benchmarking table",
  metavar = "<FILE>"
)
parser <- add_option(parser, c("-o", "--out"),
  action = "store",
  help = "Output PDF filename"
)
parser <- add_option(parser, c("-d", "--type"),
  action = "store",
  help = "Benchmarking table",
  metavar = "<>"
)
args1 <- parse_args(parser)

if (is_null(args1$table)) {
  input_dir0 <- file.path(".")
  input1 <- file.path(input_dir0, "table-some-2.tsv")
  input2 <- "memory"
  input2 <- "time"
  output1 <- file.path(input_dir0, "figure-some-2.pdf")

  args1 <- parse_args(parser, args = c("--table", input1,
                                       "-o", output1,
                                       "--type", input2))
}

output1 <- paste0(args1$out)

# Load data
df <- read_tsv(args1$table)

if (args1$type == "memory") {
# Define method order and labels
method_order <- c(
  "genomesize", "msbwt", "msbwt_polish", "nextdenovo_polish",
  "getorganelle", "ptgaul", "polap_disassemble_default", "polap_disassemble_simple",
  "oatk_nextdenovo_30", "oatk_nextdenovo_20", 
  "tippo_nextdenovo_onthq", "tippo_nextdenovo_ont",
  "pmat_nextdenovo_0.1", "pmat_nextdenovo_1.0"
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

methods_to_remove <- c(
                       "polap_disassemble_simple",
                       "polap_disassemble_check_default", 
                       "polap_disassemble_check_simple")
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
  geom_boxplot(fill = "gray80", color = "black") +  # Light gray fill
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

# Helper to convert time to minutes
convert_time_to_minutes <- function(x) {
  if (is.na(x)) return(NA_real_)
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
  "genomesize", "msbwt", "msbwt_polish", "nextdenovo_polish",
  "getorganelle", "ptgaul", "polap_disassemble_default", "polap_disassemble_simple",
  "oatk_nextdenovo_30", "oatk_nextdenovo_20", 
  "tippo_nextdenovo_onthq", "tippo_nextdenovo_ont",
  "pmat_nextdenovo_0.1", "pmat_nextdenovo_1.0"
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
methods_to_remove <- c("nextdenovo_polish",
                       "polap_disassemble_simple",
                       "polap_disassemble_check_default", 
                       "polap_disassemble_check_simple")
} else {
methods_to_remove <- c("polap_disassemble_simple",
                       "polap_disassemble_check_default", 
                       "polap_disassemble_check_simple")
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
  geom_boxplot(fill = "gray80", color = "black") +  # Light gray fill
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


