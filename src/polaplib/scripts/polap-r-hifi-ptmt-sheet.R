#!/usr/bin/env Rscript
# polap-r-hifi-ptmt-sheet.R
# Version : v1.1.8  (2025-12-02)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Draw mitogenome assembly graphs from 4 pipelines (PMAT, TIPPo, HiMT, Oatk)
# per species per row, possibly over multiple pages.
#
# Input CSV (--list) columns:
#   species, code2, pmat_png, tippo_png, himt_png, oatk_png
#

library(optparse)
suppressPackageStartupMessages({
  library(readr)
  library(png)
  library(grid)
  library(grDevices)
  library(dplyr)
})

# ------------------------------------------------------------------#
# CLI
# ------------------------------------------------------------------#
opt <- OptionParser(option_list = list(
  make_option(c("--list"),
    type = "character",
    help = "CSV with columns: species, code2, pmat_png, tippo_png, himt_png, oatk_png"
  ),
  make_option(c("--out"),
    type = "character",
    help = "Output PDF file"
  ),
  make_option(c("--title"),
    type = "character",
    default = "Mitogenome assembly graphs"
  ),
  make_option(c("--subtitle"), type = "character", default = ""),
  make_option(c("--page-width-in"), type = "double", default = 11),
  make_option(c("--page-height-in"), type = "double", default = 17),
  make_option(c("--margin-in"), type = "double", default = 0.5),
  make_option(c("--gap-mm"), type = "double", default = 4.0),
  make_option(c("--label-cex"), type = "double", default = 1.2),
  make_option(c("--title-cex"), type = "double", default = 1.0),
  make_option(c("--rows-per-page"),
    type = "integer", default = 6,
    help = "Species rows per page (default: 6)"
  ),
  make_option(c("--start-page"),
    type = "integer", default = 1,
    help = "Starting page number for footer (default: 1)"
  )
)) |>
  parse_args()

if (is.null(opt$list) || is.null(opt$out)) {
  stop("Need --list and --out")
}

# ------------------------------------------------------------------#
# Data
# ------------------------------------------------------------------#
df <- suppressMessages(read_csv(
  opt$list,
  show_col_types = FALSE,
  col_types = cols(.default = col_character())
))

required_cols <- c("species", "code2", "pmat_png", "tippo_png", "himt_png", "oatk_png")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in --list CSV: ", paste(missing_cols, collapse = ", "))
}
df <- df[, required_cols]

pipelines <- c("pmat_png", "tippo_png", "himt_png", "oatk_png")
pipeline_labels <- c("PMAT", "TIPPo", "HiMT", "Oatk")

n_species <- nrow(df)
rows_per_page <- max(1L, opt$`rows-per-page`)
n_pages <- ceiling(n_species / rows_per_page)

gapx <- unit(opt$`gap-mm`, "mm")
gapy <- unit(opt$`gap-mm`, "mm")

# ------------------------------------------------------------------#
# Page drawing
# ------------------------------------------------------------------#
draw_page <- function(df_page, page_index, n_pages) {
  n_page_rows <- nrow(df_page)
  n_cols <- length(pipelines)

  grid.newpage()

  # Outer viewport with margins
  m <- unit(opt$`margin-in`, "in")
  outer_vp <- viewport(
    x = 0.5, y = 0.5,
    width = unit(1, "npc") - 2 * m,
    height = unit(1, "npc") - 2 * m
  )
  pushViewport(outer_vp)

  # ----------------------- Title & subtitle ----------------------- #
  grid.text(
    opt$title,
    x = 0.5,
    y = unit(1, "npc"),
    just = c("center", "top"),
    gp = gpar(fontface = "bold", cex = opt$`title-cex`)
  )

  has_subtitle <- nzchar(opt$subtitle)
  if (has_subtitle) {
    grid.text(
      opt$subtitle,
      x = 0.5,
      y = unit(1, "npc") - unit(1.2, "lines"),
      just = c("center", "top"),
      gp = gpar(cex = 0.9)
    )
  }

  # Top y for column headers area
  content_top <- unit(1, "npc") - unit(if (has_subtitle) 2.4 else 1.6, "lines")

  # ----------------------- Column headers ------------------------ #
  header_vp <- viewport(
    y = content_top,
    height = unit(1.0, "lines"),
    width = unit(1, "npc"),
    just = c("center", "top")
  )
  pushViewport(header_vp)
  xs <- seq(0.1, 0.9, length.out = n_cols)
  for (i in seq_len(n_cols)) {
    grid.text(
      pipeline_labels[i],
      x = xs[i],
      y = 0,
      just = c("center", "center"),
      gp = gpar(fontface = "bold", cex = 1.0)
    )
  }
  popViewport() # header_vp

  # ----------------------- Body layout (rows) --------------------- #
  footer_height <- unit(1.2, "lines")
  header_body_gap <- unit(1.0, "lines") # small gap BETWEEN header and first row
  body_height <- content_top - footer_height - header_body_gap

  body_vp <- viewport(
    y = footer_height,
    height = body_height,
    width = unit(1, "npc"),
    just = c("center", "bottom")
  )
  pushViewport(body_vp)

  pushViewport(viewport(
    layout = grid.layout(
      nrow = rows_per_page,
      ncol = n_cols,
      widths = unit(rep(1, n_cols), "null"),
      heights = unit(rep(1, rows_per_page), "null")
    )
  ))

  for (row_i in seq_len(rows_per_page)) {
    if (row_i > n_page_rows) {
      next
    }
    sp <- df_page$species[row_i]
    cd2 <- df_page$code2[row_i]
    label_text <- sprintf("%s: %s", cd2, gsub("_", " ", sp, fixed = TRUE))

    for (col_i in seq_len(n_cols)) {
      img_path <- df_page[[pipelines[col_i]]][row_i]

      pushViewport(viewport(
        layout.pos.row = row_i,
        layout.pos.col = col_i
      ))

      pushViewport(viewport(
        x = 0.5, y = 0.5,
        width = unit(1, "npc") - 2 * gapx,
        height = unit(1, "npc") - 2 * gapy,
        clip = "off"
      ))

      # Species label ABOVE images, first column only
      if (col_i == 1L) {
        grid.text(
          label_text,
          x = 0,
          y = 0.95, # slightly below top to avoid touching header
          just = c("left", "top"),
          gp = gpar(cex = opt$`label-cex`, fontface = 3)
        )
      }

      # Image region lower in the cell so it doesn't collide with header
      img_height_frac <- 0.7
      img_y_center <- 0.35

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

      popViewport() # inner
      popViewport() # cell
    }
  }

  popViewport() # grid.layout
  popViewport() # body_vp

  # ----------------------- Footer (page number) ------------------- #
  page_no <- opt$`start-page` + page_index - 1L
  last_page <- opt$`start-page` + n_pages - 1L
  label <- sprintf("page %d/%d", page_no, last_page)

  grid.text(
    label,
    x = 0.5,
    y = unit(0.3, "lines"),
    just = c("center", "bottom"),
    gp = gpar(cex = 0.8)
  )

  popViewport() # outer_vp
}

# ------------------------------------------------------------------#
# Generate PDF
# ------------------------------------------------------------------#
pdf(
  opt$out,
  width  = opt$`page-width-in`,
  height = opt$`page-height-in`
)

for (page_index in seq_len(n_pages)) {
  start_idx <- (page_index - 1L) * rows_per_page + 1L
  end_idx <- min(start_idx + rows_per_page - 1L, n_species)
  df_page <- df[start_idx:end_idx, , drop = FALSE]
  draw_page(df_page, page_index, n_pages)
}

dev.off()
cat(opt$out, "\n")
