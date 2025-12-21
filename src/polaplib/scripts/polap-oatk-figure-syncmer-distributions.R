#!/usr/bin/env Rscript
###############################################################################
# scripts/polap-oatk-figure-syncmer-distributions.R
#
# Version: v0.1.2
# Author : Sang Chul Choi (POLAP)
# License: GPL-3.0+
#
# Purpose:
#   Draw conceptual schematic figures for:
#     (A) Subsampling effect
#     (B) Tiara filtering effect
#
#   Outputs:
#     1) Annotated versions (for standalone viewing):
#        <out_prefix>-A.png
#        <out_prefix>-B.png
#     2) Clean versions (for embedding into a combined figure):
#        <out_prefix>-A.clean.png
#        <out_prefix>-B.clean.png
#
# Usage:
#   Rscript --vanilla scripts/polap-oatk-figure-syncmer-distributions.R \
#     --out-prefix oatk-syncmer-schematic \
#     --dpi 600 --width-in 10 --height-in 4
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(optparse)
})

opt_list <- list(
  make_option(c("-o", "--out-prefix"),
    type = "character", default = "oatk-syncmer-schematic",
    help = "Output prefix (default: %default)"
  ),
  make_option(c("--dpi"),
    type = "integer", default = 600,
    help = "PNG DPI (default: %default)"
  ),
  make_option(c("--width-in"),
    type = "double", default = 10,
    help = "Width in inches (default: %default)"
  ),
  make_option(c("--height-in"),
    type = "double", default = 4,
    help = "Height in inches (default: %default)"
  ),
  make_option(c("--c-moderate"),
    type = "double", default = 10,
    help = "Moderate cutoff line position (default: %default)"
  ),
  make_option(c("--c-low"),
    type = "double", default = 5,
    help = "Low cutoff line position (default: %default)"
  ),
  make_option(c("--peak"),
    type = "double", default = 15,
    help = "Peak location of conceptual nuclear distribution (default: %default)"
  ),
  make_option(c("--sd"),
    type = "double", default = 6,
    help = "SD of conceptual nuclear distribution (default: %default)"
  ),
  make_option(c("--tiara-scale"),
    type = "double", default = 0.4,
    help = "Scale factor for 'After Tiara' distribution (default: %default)"
  ),
  make_option(c("--title-size"),
    type = "double", default = 18,
    help = "Title font size for clean figures (default: %default)"
  )
)

opt <- parse_args(OptionParser(option_list = opt_list))

# ---- helper: save high-quality PNG ----
# Prefer ragg (crisp text), fallback to base png.
save_png <- function(filename, plot, width_in, height_in, dpi) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename,
      width = width_in, height = height_in,
      units = "in", res = dpi, scaling = 1
    )
    print(plot)
    dev.off()
  } else {
    grDevices::png(filename,
      width = width_in, height = height_in,
      units = "in", res = dpi
    )
    print(plot)
    dev.off()
  }
}

# ---- data: conceptual distribution ----
x <- seq(0, 50, length.out = 501)
y_nuc <- dnorm(x, mean = opt$peak, sd = opt$sd)
y_nuc <- y_nuc / max(y_nuc) # normalize peak to 1

df <- data.frame(
  x = x,
  y_nuc = y_nuc,
  y_tiara = opt$`tiara-scale` * y_nuc
)

c_mod <- opt$`c-moderate`
c_low <- opt$`c-low`

# Shared for Panel A
dfA_all <- data.frame(x = df$x, y = df$y_nuc)
dfA_removed_mod <- subset(dfA_all, x < c_mod)

# ------------------------------------------------------------------
# Panel A (annotated; unchanged behavior)
# ------------------------------------------------------------------
pA <- ggplot() +
  geom_area(data = dfA_all, aes(x = x, y = y), alpha = 0.90, fill = "grey80") +
  geom_area(data = dfA_removed_mod, aes(x = x, y = y), alpha = 0.60, fill = "tomato") +
  geom_vline(xintercept = c_mod, linetype = "dashed", linewidth = 0.7) +
  annotate("text",
    x = c_mod, y = 1.05, label = "c ≈ 8–11",
    vjust = 0, hjust = 0.5, size = 4
  ) +
  annotate("text",
    x = 3.5, y = 0.35,
    label = "Removed low-coverage\nnuclear background",
    hjust = 0, vjust = 0.5, size = 4
  ) +
  annotate("text",
    x = 22, y = 0.62,
    label = "Kept syncmers\n(≥ c)",
    hjust = 0, vjust = 0.5, size = 4
  ) +
  labs(
    title = "(A) Subsampling effect",
    x = "Syncmer coverage",
    y = "Relative abundance"
  ) +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 1.2), expand = FALSE) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------
# Panel B (annotated; with left-of-cutoff coloring at c = 5)
# ------------------------------------------------------------------
dfB_unf <- dfA_all
dfB_tia <- data.frame(x = df$x, y = df$y_tiara)

dfB_unf_left <- subset(dfB_unf, x < c_low)
dfB_tia_left <- subset(dfB_tia, x < c_low)

pB <- ggplot() +
  # full unfiltered baseline
  geom_area(data = dfB_unf, aes(x = x, y = y), fill = "grey80", alpha = 0.90) +
  # full Tiara overlay
  geom_area(data = dfB_tia, aes(x = x, y = y), fill = "#f8766d", alpha = 0.40) +
  # highlight left-of-c=5
  geom_area(data = dfB_unf_left, aes(x = x, y = y), fill = "grey50", alpha = 0.60) +
  geom_area(data = dfB_tia_left, aes(x = x, y = y), fill = "#f8766d", alpha = 0.80) +
  geom_vline(xintercept = c_low, linetype = "dotted", linewidth = 0.7) +
  geom_vline(xintercept = c_mod, linetype = "dashed", linewidth = 0.7) +
  annotate("text",
    x = c_low, y = 1.05, label = "c = 5",
    vjust = 0, hjust = 0.5, size = 4
  ) +
  annotate("text",
    x = c_mod, y = 1.05, label = "c ≈ 8–11",
    vjust = 0, hjust = 0.5, size = 4
  ) +
  annotate("label",
    x = 28, y = 0.95,
    label = "High background\nwithout filtering",
    hjust = 0, vjust = 1, size = 4,
    label.size = 0.2, alpha = 0.90
  ) +
  annotate("label",
    x = 28, y = 0.20,
    label = "Lower cutoff usable\nwith Tiara",
    hjust = 0, vjust = 0, size = 4,
    label.size = 0.2, alpha = 0.90
  ) +
  labs(
    title = "(B) Tiara filtering effect",
    x = "Syncmer coverage",
    y = "Relative abundance"
  ) +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 1.2), expand = FALSE) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  )

# ------------------------------------------------------------------
# CLEAN versions (axes + cutoff labels; left-aligned title)
# ------------------------------------------------------------------

clean_theme <- function(title_size) {
  theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0, # LEFT align title
        size = title_size
      ),
      axis.title.x = element_text(size = 24),
      axis.title.y = element_text(size = 24),
      panel.grid.minor = element_blank()
    )
}

dfA_left_low <- subset(dfA_all, x < c_low)
dfA_left_mod <- subset(dfA_all, x < c_mod)

# ---------- Panel A (clean) ----------
pA_clean <- ggplot() +
  # full baseline
  geom_area(
    data = dfA_all, aes(x = x, y = y),
    alpha = 0.90, fill = "grey80"
  ) +
  # left of c = 10 (moderate cutoff)
  geom_area(
    data = dfA_left_mod, aes(x = x, y = y),
    alpha = 0.50, fill = "tomato"
  ) +
  # left of c = 5 (low cutoff, darker)
  geom_area(
    data = dfA_left_low, aes(x = x, y = y),
    alpha = 0.80, fill = "tomato"
  ) +
  geom_vline(xintercept = c_low, linetype = "dotted", linewidth = 0.7) +
  geom_vline(xintercept = c_mod, linetype = "dashed", linewidth = 0.7) +
  annotate("text",
    x = c_low, y = 1.05,
    label = "c = 5", vjust = 0, hjust = 0.5, size = 4
  ) +
  annotate("text",
    x = c_mod, y = 1.05,
    label = "c = 10", vjust = 0, hjust = 0.5, size = 4
  ) +
  labs(
    title = "(A) Subsampling effect",
    x = "Syncmer coverage",
    y = "Relative abundance"
  ) +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 1.2), expand = FALSE) +
  clean_theme(opt$`title-size`)

# ---------- Panel B (clean) ----------
pB_clean <- ggplot() +
  geom_area(
    data = dfB_unf, aes(x = x, y = y),
    fill = "grey80", alpha = 0.90
  ) +
  geom_area(
    data = dfB_tia, aes(x = x, y = y),
    fill = "#f8766d", alpha = 0.40
  ) +
  geom_area(
    data = dfB_unf_left, aes(x = x, y = y),
    fill = "grey50", alpha = 0.60
  ) +
  geom_area(
    data = dfB_tia_left, aes(x = x, y = y),
    fill = "#f8766d", alpha = 0.80
  ) +
  geom_vline(xintercept = c_low, linetype = "dotted", linewidth = 0.7) +
  geom_vline(xintercept = c_mod, linetype = "dashed", linewidth = 0.7) +
  annotate("text",
    x = c_low, y = 1.05,
    label = "c = 5", vjust = 0, hjust = 0.5, size = 4
  ) +
  annotate("text",
    x = c_mod, y = 1.05,
    label = "c = 10", vjust = 0, hjust = 0.5, size = 4
  ) +
  labs(
    title = "(B) Prefiltering effect",
    x = "Syncmer coverage",
    y = "Relative abundance"
  ) +
  coord_cartesian(xlim = c(0, 50), ylim = c(0, 1.2), expand = FALSE) +
  clean_theme(opt$`title-size`)


# ---- Output ----
outA <- sprintf("%s-A.png", opt$`out-prefix`)
outB <- sprintf("%s-B.png", opt$`out-prefix`)
outA_clean <- sprintf("%s-A.clean.png", opt$`out-prefix`)
outB_clean <- sprintf("%s-B.clean.png", opt$`out-prefix`)

save_png(outA, pA, opt$`width-in`, opt$`height-in`, opt$dpi)
save_png(outB, pB, opt$`width-in`, opt$`height-in`, opt$dpi)
save_png(outA_clean, pA_clean, opt$`width-in`, opt$`height-in`, opt$dpi)
save_png(outB_clean, pB_clean, opt$`width-in`, opt$`height-in`, opt$dpi)

message(sprintf("[OK] wrote %s", outA))
message(sprintf("[OK] wrote %s", outB))
message(sprintf("[OK] wrote %s", outA_clean))
message(sprintf("[OK] wrote %s", outB_clean))
