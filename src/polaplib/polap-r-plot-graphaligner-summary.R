#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)

summary_tsv <- args[1] # e.g., summary.tsv
output_pdf <- args[2] # e.g., summary.pdf

library(ggplot2)
library(readr)
library(dplyr)

df <- read_tsv(summary_tsv)

df$keep <- factor(df$keep, levels = c("yes", "no"), labels = c("Retained", "Removed"))

p1 <- ggplot(df, aes(x = identity, fill = keep)) +
  geom_density(alpha = 0.6) +
  labs(title = "Identity distribution", x = "Identity", y = "Density")

p2 <- ggplot(df, aes(x = left_clip, fill = keep)) +
  geom_density(alpha = 0.6) +
  labs(title = "Left clip distribution", x = "Left soft clip (bp)", y = "Density")

p3 <- ggplot(df, aes(x = right_clip, fill = keep)) +
  geom_density(alpha = 0.6) +
  labs(title = "Right clip distribution", x = "Right soft clip (bp)", y = "Density")

pdf(output_pdf, width = 9, height = 8)
print(p1)
print(p2)
print(p3)
dev.off()
