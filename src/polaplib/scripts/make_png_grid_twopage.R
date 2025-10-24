#!/usr/bin/env Rscript
# make_png_grid_twopage.R
# Version : v1.2.0  (2025-10-07)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Render a two-page PDF: page 1 from PT list, page 2 from MT list.
# Each list is a CSV with columns: species,png (png may be empty).
# Missing PNG -> bold red "No assembly" on white background.
# Uses same typographic/gap controls as make_png_grid.R v1.1.0.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(png)
  library(grid)
  library(grDevices)
})

opt <- OptionParser(option_list = list(
  make_option(c("--pt-list"), type = "character", help = "CSV species,png for PT"),
  make_option(c("--mt-list"), type = "character", help = "CSV species,png for MT"),
  make_option(c("--out"), type = "character", help = "Output PDF"),
  make_option(c("--pt-title"), type = "character", default = "Plastid assemblies (pt)"),
  make_option(c("--mt-title"), type = "character", default = "Mitochondrial assemblies (mt)"),
  make_option(c("--page-width-in"), type = "double", default = 11),
  make_option(c("--page-height-in"), type = "double", default = 8.5),
  make_option(c("--margin-in"), type = "double", default = 0.45),
  make_option(c("--rows"), type = "integer", default = 0),
  make_option(c("--cols"), type = "integer", default = 0),
  make_option(c("--gap-mm"), type = "double", default = 4.0),
  make_option(c("--label-cex"), type = "double", default = 0.60),
  make_option(c("--title-cex"), type = "double", default = 1.0),
  make_option(c("--label-strip-frac"), type = "double", default = 0.16)
)) |> parse_args()

stopifnot(!is.null(opt$`pt-list`), !is.null(opt$`mt-list`), !is.null(opt$out))

read_list <- function(path) {
  d <- suppressMessages(read_csv(path,
    show_col_types = FALSE,
    col_types = cols(
      species = col_character(),
      png = col_character()
    )
  ))
  if (!all(c("species", "png") %in% names(d))) stop("[ERR] CSV must have species,png columns: ", path)
  d
}

render_page <- function(df, title_text) {
  N <- nrow(df)
  if (N == 0) df <- data.frame(species = character(), png = character())
  rows <- as.integer(opt$rows)
  cols <- as.integer(opt$cols)
  if (rows <= 0 && cols <= 0) {
    cols <- ceiling(sqrt(N))
    rows <- ceiling(ifelse(N == 0, 1, ceiling(N / cols)))
  } else if (rows <= 0 && cols > 0) {
    rows <- ceiling(ifelse(N == 0, 1, ceiling(N / cols)))
  } else if (rows > 0 && cols <= 0) {
    cols <- ceiling(ifelse(N == 0, 1, ceiling(N / rows)))
  } else if (rows * cols < N) {
    rows <- ceiling(ifelse(cols == 0, 1, ceiling(N / max(1, cols))))
  }

  grid.newpage()

  m <- unit(opt$`margin-in`, "in")
  pushViewport(viewport(
    x = .5, y = .5,
    width = unit(1, "npc") - 2 * m,
    height = unit(1, "npc") - 2 * m
  ))
  grid.text(title_text,
    x = .5, y = unit(1, "npc"), just = c("center", "top"),
    gp = gpar(cex = opt$`title-cex`, fontface = "bold")
  )

  pushViewport(viewport(
    y = unit(0, "npc"), just = "bottom",
    width = unit(1, "npc"), height = unit(0.88, "npc")
  ))

  gapx <- unit(opt$`gap-mm`, "mm")
  gapy <- unit(opt$`gap-mm`, "mm")
  pushViewport(viewport(layout = grid.layout(
    nrow = rows, ncol = cols,
    widths = unit(rep(1, cols), "null"),
    heights = unit(rep(1, rows), "null")
  )))

  lbl_frac <- max(0.10, min(opt$`label-strip-frac`, 0.28))
  img_height_frac <- 1 - lbl_frac - 0.04
  img_y_center <- lbl_frac + img_height_frac / 2 + 0.02

  k <- 1
  for (r in seq_len(rows)) {
    for (c in seq_len(cols)) {
      if (k > N) break
      sp <- df$species[k]
      img_path <- df$png[k]

      pushViewport(viewport(layout.pos.row = r, layout.pos.col = c))
      pushViewport(viewport(
        x = 0.5, y = 0.5,
        width = unit(1, "npc") - 2 * gapx,
        height = unit(1, "npc") - 2 * gapy,
        clip = "on"
      ))

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
}

# main
pt <- read_list(opt$`pt-list`)
mt <- read_list(opt$`mt-list`)

pdf(opt$out, width = opt$`page-width-in`, height = opt$`page-height-in`)
render_page(pt, opt$`pt-title`)
render_page(mt, opt$`mt-title`)
dev.off()
cat(opt$out, "\n")
