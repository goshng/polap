#!/usr/bin/env Rscript
# scripts/polap-r-hifi-ptmt-sheet.R
# Version : v1.1.0  (2025-12-02)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Draw mitogenome assembly graphs from 4 pipelines (pmat, tippo, himt, oatk)
# per species per row.
#
# Input CSV (via --list) must have columns:
#   species, pmat_png, tippo_png, himt_png, oatk_png
#
# For each species:
#   row = species
#   columns = [PMAT, TIPPo, HiMT, Oatk]
# If a PNG does not exist or is empty, "No assembly" is shown.

library(optparse)
suppressPackageStartupMessages({
  library(readr)
  library(png)
  library(grid)
  library(grDevices)
  library(dplyr)
})

opt <- OptionParser(option_list = list(
  make_option(c("--list"),
    type = "character",
    help = "CSV with columns: species, pmat_png, tippo_png, himt_png, oatk_png"
  ),
  make_option(c("--out"),
    type = "character",
    help = "Output PDF file"
  ),
  make_option(c("--title"),
    type = "character",
    default = "Mitogenome assembly graphs (PMAT, TIPPo, HiMT, Oatk)"
  ),
  make_option(c("--page-width-in"), type = "double", default = 11),
  make_option(c("--page-height-in"), type = "double", default = 8.5),
  make_option(c("--margin-in"), type = "double", default = 0.45),
  make_option(c("--gap-mm"), type = "double", default = 4.0),
  make_option(c("--label-cex"), type = "double", default = 0.55),
  make_option(c("--title-cex"), type = "double", default = 1.0)
)) |> parse_args()

stopifnot(!is.null(opt$list), file.exists(opt$list), !is.null(opt$out))

read_list <- function(path) {
  suppressMessages(read_csv(
    path,
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  ))
}

df <- read_list(opt$list)

# basic sanity
required_cols <- c("species", "pmat_png", "tippo_png", "himt_png", "oatk_png")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in --list CSV: ", paste(missing_cols, collapse = ", "))
}

# Keep only the required columns, preserve order
df <- df[, required_cols]

pipelines <- c("pmat_png", "tippo_png", "himt_png", "oatk_png")
pipeline_labels <- c("PMAT", "TIPPo", "HiMT", "Oatk")

draw_page <- function(df) {
  n_species <- nrow(df)
  n_pipelines <- length(pipelines)

  # Layout:
  #   row 1: column headers (pipelines)
  #   rows 2..(n_species+1): species
  n_rows_layout <- n_species + 1L
  n_cols_layout <- n_pipelines

  grid.newpage()
  m <- unit(opt$`margin-in`, "in")

  # Outer viewport with margins
  pushViewport(viewport(
    x = .5, y = .5,
    width = unit(1, "npc") - 2 * m,
    height = unit(1, "npc") - 2 * m
  ))

  # Title at top
  grid.text(
    opt$title,
    x = .5, y = unit(1, "npc"),
    just = c("center", "top"),
    gp = gpar(cex = opt$`title-cex`, fontface = "bold")
  )

  # Body viewport (slightly below title)
  pushViewport(viewport(
    y = unit(0, "npc"), just = "bottom",
    width = unit(1, "npc"),
    height = unit(0.90, "npc")
  ))

  gapx <- unit(opt$`gap-mm`, "mm")
  gapy <- unit(opt$`gap-mm`, "mm")

  pushViewport(viewport(layout = grid.layout(
    nrow    = n_rows_layout,
    ncol    = n_cols_layout,
    widths  = unit(rep(1, n_cols_layout), "null"),
    heights = unit(rep(1, n_rows_layout), "null")
  )))

  # ------------------------------------------------------------------------#
  # Header row: pipeline labels
  # ------------------------------------------------------------------------#
  for (c in seq_len(n_pipelines)) {
    pushViewport(viewport(
      layout.pos.row = 1L,
      layout.pos.col = c
    ))
    grid.text(
      pipeline_labels[c],
      x = 0.5, y = 0.95,
      just = c("center", "top"),
      gp = gpar(fontface = "bold")
    )
    popViewport()
  }

  # ------------------------------------------------------------------------#
  # Species rows
  # ------------------------------------------------------------------------#
  for (i in seq_len(n_species)) {
    sp <- df$species[i]

    for (c in seq_len(n_pipelines)) {
      col_name <- pipelines[c]
      img_path <- df[[col_name]][i]

      pushViewport(viewport(
        layout.pos.row = i + 1L, # +1 because row 1 is header
        layout.pos.col = c
      ))

      # inner viewport with gap margins
      pushViewport(viewport(
        x = .5, y = .5,
        width = unit(1, "npc") - 2 * gapx,
        height = unit(1, "npc") - 2 * gapy,
        clip = "on"
      ))

      # main image area
      img_height_frac <- 0.80
      img_y_center <- 0.50

      if (!is.na(img_path) && nzchar(img_path) && file.exists(img_path)) {
        im <- try(readPNG(img_path), silent = TRUE)
        if (!inherits(im, "try-error")) {
          grid.raster(
            im,
            width  = unit(0.95, "npc"),
            height = unit(img_height_frac, "npc"),
            y      = unit(img_y_center, "npc")
          )
        } else {
          grid.text(
            "No assembly",
            gp = gpar(col = "red3", cex = 0.9, fontface = "bold"),
            y  = unit(img_y_center, "npc")
          )
        }
      } else {
        grid.text(
          "No assembly",
          gp = gpar(col = "red3", cex = 0.9, fontface = "bold"),
          y  = unit(img_y_center, "npc")
        )
      }

      # Species label at bottom of the row, but only in the first column to avoid repetition
      if (c == 1L) {
        grid.text(
          gsub("_", " ", sp, fixed = TRUE),
          x = 0.02, y = 0.02,
          just = c("left", "bottom"),
          gp = gpar(cex = opt$`label-cex`, fontface = 3)
        )
      }

      popViewport() # inner
      popViewport() # cell
    }
  }

  popViewport() # layout
  popViewport() # body
  popViewport() # outer
}

pdf(opt$out, width = opt$`page-width-in`, height = opt$`page-height-in`)
draw_page(df)
dev.off()
cat(opt$out, "\n")
