#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(jsonlite)
  library(dplyr)
  library(ggplot2)
  suppressWarnings(requireNamespace("ggrepel", quietly = TRUE))
})

option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Path to manifest JSON"),
  make_option(c("-o", "--output"), type = "character", help = "Output figure file (PDF/PNG)"),
  make_option(c("--width"), type = "numeric", default = 8),
  make_option(c("--height"), type = "numeric", default = 6),
  make_option(c("--title"), type = "character", default = "Mitochondrial genome assemblies (v6)")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$input), !is.null(opt$output))

cat("[INFO] Reading manifest:", opt$input, "\n")
# IMPORTANT: keep nested lists; do NOT flatten
manifest <- fromJSON(opt$input, simplifyVector = FALSE)

items <- manifest$items
if (is.null(items) || length(items) == 0L) {
  stop("[ERR] No items in manifest; cannot make assembly figure.")
}

# helper to safely pull nested values
get_num <- function(x) {
  if (is.null(x)) {
    return(NA_real_)
  }
  as.numeric(x)
}

rows <- lapply(items, function(it) {
  # Only plot mt for now (extend to pt later as needed)
  mt <- it$mt
  st <- if (!is.null(mt)) mt$stats else NULL
  if (is.null(st)) {
    return(NULL)
  }

  data.frame(
    species = as.character(it$species %||% NA_character_),
    organelle = "mt",
    total_len = get_num(st$total_len),
    N50 = get_num(st$N50),
    n_segments = get_num(st$n_segments),
    max_seg = get_num(st$max_seg),
    is_circular = factor(ifelse(get_num(st$is_circular) == 1, "circular", "linear"),
      levels = c("linear", "circular")
    ),
    check.names = FALSE
  )
})

rows <- rows[!vapply(rows, is.null, logical(1))]
if (length(rows) == 0L) stop("[ERR] No mt.stats blocks found in manifest.")

df <- bind_rows(rows)

# order species by genome size descending
df$species <- factor(df$species, levels = df$species[order(df$total_len, decreasing = TRUE)])

# build plot
p <- ggplot(df, aes(x = species, y = total_len / 1e3, fill = is_circular)) +
  geom_col() +
  {
    if ("ggrepel" %in% loadedNamespaces()) {
      ggrepel::geom_text_repel(aes(label = sprintf("N50=%.0fk", N50 / 1e3)),
        size = 3, color = "black", direction = "y", max.overlaps = Inf
      )
    } else {
      geom_text(aes(label = sprintf("N50=%.0fk", N50 / 1e3)),
        size = 3, color = "black", angle = 90, vjust = -0.2
      )
    }
  } +
  labs(
    title = opt$title,
    x = "Species",
    y = "Total length (kb)",
    fill = "Topology"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title  = element_text(hjust = 0.5)
  )

ext <- tolower(tools::file_ext(opt$output))
if (ext == "pdf") {
  ggsave(opt$output, plot = p, width = opt$width, height = opt$height, device = "pdf")
} else {
  ggsave(opt$output, plot = p, width = opt$width, height = opt$height, dpi = 300)
}
cat("[OK] Wrote figure:", opt$output, "\n")
