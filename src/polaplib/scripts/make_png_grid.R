#!/usr/bin/env Rscript
# make_png_grid.R
# Version : v1.1.0  (2025-10-07)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Single-page grid of assembly PNGs with true inter-panel gaps.
# Missing PNGs show bold red "No assembly" on white background.
# Typography is configurable: label size, title size, optional no-title.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(png)
  library(grid)
  library(grDevices)
})

opt <- OptionParser(option_list = list(
  make_option(c("--list"), type = "character", help = "CSV with columns: species,png"),
  make_option(c("--out"), type = "character", help = "Output PDF"),
  make_option(c("--title"), type = "character", default = "Assemblies"),
  make_option(c("--no-title"),
    action = "store_true", default = FALSE,
    help = "Suppress the main title"
  ),
  make_option(c("--title-cex"),
    type = "double", default = 1.0,
    help = "Title size multiplier [default 1.0]"
  ),
  make_option(c("--label-cex"),
    type = "double", default = 0.60,
    help = "Species label size multiplier [default 0.60]"
  ),
  make_option(c("--page-width-in"), type = "double", default = 8.5),
  make_option(c("--page-height-in"), type = "double", default = 11),
  make_option(c("--margin-in"), type = "double", default = 0.45),
  make_option(c("--rows"), type = "integer", default = 0),
  make_option(c("--cols"), type = "integer", default = 0),
  make_option(c("--gap-mm"),
    type = "double", default = 4.0,
    help = "Gap between panels in millimetres [default 4.0]"
  ),
  make_option(c("--label-strip-frac"),
    type = "double", default = 0.16,
    help = "Fractional height for label strip [default 0.16]"
  )
)) |> parse_args()

stopifnot(!is.null(opt$list), !is.null(opt$out))
df <- suppressMessages(
  read_csv(opt$list,
    show_col_types = FALSE,
    col_types = cols(species = col_character(), png = col_character())
  )
)
stopifnot(all(c("species", "png") %in% names(df)))
if (nrow(df) == 0) stop("[ERR] empty list")

N <- nrow(df)
rows <- as.integer(opt$rows)
cols <- as.integer(opt$cols)

# single-line else-if (R syntax)
if (rows <= 0 && cols <= 0) {
  cols <- ceiling(sqrt(N))
  rows <- ceiling(N / cols)
} else if (rows <= 0 && cols > 0) {
  rows <- ceiling(N / cols)
} else if (rows > 0 && cols <= 0) {
  cols <- ceiling(N / rows)
} else if (rows * cols < N) {
  rows <- ceiling(N / cols)
}

# PDF device
pdf(opt$out, width = opt$`page-width-in`, height = opt$`page-height-in`)
grid.newpage()

m <- unit(opt$`margin-in`, "in")
pushViewport(viewport(
  x = .5, y = .5,
  width = unit(1, "npc") - 2 * m,
  height = unit(1, "npc") - 2 * m
))

# Title (optional)
if (!isTRUE(opt$`no-title`)) {
  grid.text(opt$title,
    x = .5, y = unit(1, "npc"), just = c("center", "top"),
    gp = gpar(cex = opt$`title-cex`, fontface = "bold")
  )
}

# Content area: keep ~12% for title line if shown; else use full height
content_height_frac <- if (isTRUE(opt$`no-title`)) 1.0 else 0.88
pushViewport(viewport(
  y = unit(0, "npc"), just = "bottom",
  width = unit(1, "npc"), height = unit(content_height_frac, "npc")
))

# Real gaps
gapx <- unit(opt$`gap-mm`, "mm")
gapy <- unit(opt$`gap-mm`, "mm")

pushViewport(viewport(layout = grid.layout(
  nrow = rows, ncol = cols,
  widths = unit(rep(1, cols), "null"),
  heights = unit(rep(1, rows), "null")
)))

# Panel layout constants
lbl_frac <- max(0.10, min(opt$`label-strip-frac`, 0.28)) # a bit taller than before
img_height_frac <- 1 - lbl_frac - 0.04 # image frac + cushion
img_y_center <- lbl_frac + img_height_frac / 2 + 0.02 # lift image above label

k <- 1
for (r in seq_len(rows)) {
  for (c in seq_len(cols)) {
    if (k > N) break
    sp <- df$species[k]
    img_path <- df$png[k]

    # Cell -> inner padded viewport = visible gap
    pushViewport(viewport(layout.pos.row = r, layout.pos.col = c))
    pushViewport(viewport(
      x = 0.5, y = 0.5,
      width = unit(1, "npc") - 2 * gapx,
      height = unit(1, "npc") - 2 * gapy,
      clip = "on"
    ))

    # Image or "No assembly"
    if (!is.na(img_path) && nzchar(img_path) && file.exists(img_path)) {
      img <- try(readPNG(img_path), silent = TRUE)
      if (!inherits(img, "try-error")) {
        grid.raster(img,
          width  = unit(0.95, "npc"),
          height = unit(img_height_frac, "npc"),
          y      = unit(img_y_center, "npc")
        )
      } else {
        grid.text("No assembly",
          gp = gpar(col = "red3", cex = .95, fontface = "bold"),
          y = unit(img_y_center, "npc")
        )
      }
    } else {
      grid.text("No assembly",
        gp = gpar(col = "red3", cex = .95, fontface = "bold"),
        y = unit(img_y_center, "npc")
      )
    }

    # Label strip (taller) + species name smaller and clearly separated
    grid.rect(
      x = unit(0.5, "npc"),
      y = unit(lbl_frac / 2, "npc"),
      width = unit(0.98, "npc"),
      height = unit(lbl_frac, "npc"),
      gp = gpar(fill = "white", col = NA)
    )
    grid.text(gsub("_", " ", sp, fixed = TRUE),
      x = unit(0.5, "npc"), y = unit(0.00, "npc"),
      just = c("center", "bottom"),
      gp = gpar(cex = opt$`label-cex`, fontface = 3)
    )

    popViewport()
    popViewport()
    k <- k + 1
  }
}
popViewport()
popViewport()
dev.off()
cat(opt$out, "\n")
