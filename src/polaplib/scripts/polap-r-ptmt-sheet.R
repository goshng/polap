#!/usr/bin/env Rscript
# polap-r-ptmt-sheet.R
# Version : v1.0.1  (2025-10-07)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+

library(optparse)
suppressPackageStartupMessages({
  library(readr)
  library(png)
  library(grid)
  library(grDevices)
  library(dplyr)
})

`%||%` <- function(a, b) if (is.null(a) || is.na(a) || !nzchar(a)) b else a

opt <- OptionParser(option_list = list(
  make_option(c("--pt-list"), type = "character"),
  make_option(c("--mt-list"), type = "character"),
  make_option(c("--annot-pt"), type = "character"),
  make_option(c("--annot-mt"), type = "character"),
  make_option(c("--out"), type = "character"),
  make_option(c("--pt-title"), type = "character", default = "Plastid assemblies (pt)"),
  make_option(c("--mt-title"), type = "character", default = "Mitochondrial assemblies (mt)"),
  make_option(c("--page-width-in"), type = "double", default = 11),
  make_option(c("--page-height-in"), type = "double", default = 8.5),
  make_option(c("--margin-in"), type = "double", default = 0.45),
  make_option(c("--rows"), type = "integer", default = 0),
  make_option(c("--cols"), type = "integer", default = 0),
  make_option(c("--gap-mm"), type = "double", default = 4.0),
  make_option(c("--label-cex"), type = "double", default = 0.55),
  make_option(c("--title-cex"), type = "double", default = 1.0),
  make_option(c("--label-strip-frac"), type = "double", default = 0.16)
)) |> parse_args()

stopifnot(!is.null(opt$`pt-list`), !is.null(opt$`mt-list`), !is.null(opt$out))

read_list <- function(path) {
  suppressMessages(read_csv(path,
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  ))
}

# Read an annotation CSV and **normalize columns** to:
#   species, mt_genes, pt_genes, genome_kb   (all character; we’ll coerce later)
read_annot <- function(path, org = c("pt", "mt")) {
  org <- match.arg(org)
  if (is.null(path) || !file.exists(path)) {
    return(NULL)
  }
  d <- suppressMessages(read_csv(path,
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  ))
  stopifnot("species" %in% names(d))
  # If legacy schema {species, genes, genome_kb}, map genes → *_{org}_genes
  if (!("mt_genes" %in% names(d)) && !("pt_genes" %in% names(d)) && ("genes" %in% names(d))) {
    if (org == "mt") d$mt_genes <- d$genes else d$pt_genes <- d$genes
  }
  if (!("mt_genes" %in% names(d))) d$mt_genes <- NA_character_
  if (!("pt_genes" %in% names(d))) d$pt_genes <- NA_character_
  if (!("genome_kb" %in% names(d))) d$genome_kb <- NA_character_
  d
}

# Draw one page (org = "pt" or "mt"); overlay appropriate gene totals
draw_page <- function(df, annot, title_text, org = c("pt", "mt")) {
  org <- match.arg(org)
  N <- nrow(df)
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
  pushViewport(viewport(x = .5, y = .5, width = unit(1, "npc") - 2 * m, height = unit(1, "npc") - 2 * m))
  grid.text(title_text,
    x = .5, y = unit(1, "npc"), just = c("center", "top"),
    gp = gpar(cex = opt$`title-cex`, fontface = "bold")
  )

  pushViewport(viewport(y = unit(0, "npc"), just = "bottom", width = unit(1, "npc"), height = unit(0.88, "npc")))
  gapx <- unit(opt$`gap-mm`, "mm")
  gapy <- unit(opt$`gap-mm`, "mm")
  pushViewport(viewport(layout = grid.layout(
    nrow = rows, ncol = cols,
    widths = unit(rep(1, cols), "null"), heights = unit(rep(1, rows), "null")
  )))

  lbl_frac <- max(0.10, min(opt$`label-strip-frac`, 0.28))
  img_height_frac <- 1 - lbl_frac - 0.04
  img_y_center <- lbl_frac + img_height_frac / 2 + 0.02

  # join annotation
  if (!is.null(annot) && nrow(annot) > 0) {
    df <- df %>% left_join(annot, by = "species")
  } else {
    df$mt_genes <- NA_character_
    df$pt_genes <- NA_character_
    df$genome_kb <- NA_character_
  }

  # coerce numeric for formatting
  to_num <- function(x) suppressWarnings(as.numeric(x))
  df$mt_genes_n <- to_num(df$mt_genes)
  df$pt_genes_n <- to_num(df$pt_genes)
  df$genome_kb_n <- to_num(df$genome_kb)

  k <- 1
  for (r in seq_len(rows)) {
    for (c in seq_len(cols)) {
      if (k > nrow(df)) break
      sp <- df$species[k]
      img <- df$png[k]

      pushViewport(viewport(layout.pos.row = r, layout.pos.col = c))
      pushViewport(viewport(x = .5, y = .5, width = unit(1, "npc") - 2 * gapx, height = unit(1, "npc") - 2 * gapy, clip = "on"))

      if (!is.na(img) && nzchar(img) && file.exists(img)) {
        im <- try(readPNG(img), silent = TRUE)
        if (!inherits(im, "try-error")) {
          grid.raster(im, width = unit(0.95, "npc"), height = unit(img_height_frac, "npc"), y = unit(img_y_center, "npc"))
        } else {
          grid.text("No assembly", gp = gpar(col = "red3", cex = .95, fontface = "bold"), y = unit(img_y_center, "npc"))
        }
      } else {
        grid.text("No assembly", gp = gpar(col = "red3", cex = .95, fontface = "bold"), y = unit(img_y_center, "npc"))
      }

      # overlay string — choose fields by page/org
      mtg <- df$mt_genes_n[k]
      ptg <- df$pt_genes_n[k]
      gkb <- df$genome_kb_n[k]
      ann <- NA_character_
      if (org == "pt") {
        if (!is.na(gkb)) {
          ann <- sprintf("PT=%.0f, MT=%.0f, len=%.0f kb", ptg %||% NA, mtg %||% NA, gkb)
        } else {
          ann <- sprintf("PT=%.0f, MT=%.0f", ptg %||% NA, mtg %||% NA)
        }
      } else {
        if (!is.na(gkb)) {
          ann <- sprintf("len=%.0f kb, MT=%.0f, PT=%.0f", gkb, mtg %||% NA, ptg %||% NA)
        } else {
          ann <- sprintf("MT=%.0f, PT=%.0f", mtg %||% NA, ptg %||% NA)
        }
      }
      if (!is.na(ann)) {
        grid.text(ann, x = unit(0.02, "npc"), y = unit(0.98, "npc"), just = c("left", "top"), gp = gpar(col = "black", cex = 0.55))
      }

      grid.rect(
        x = unit(0.5, "npc"), y = unit(lbl_frac / 2, "npc"),
        width = unit(0.98, "npc"), height = unit(lbl_frac, "npc"),
        gp = gpar(fill = "white", col = NA)
      )
      grid.text(gsub("_", " ", sp, fixed = TRUE),
        x = unit(0.5, "npc"), y = unit(0.00, "npc"),
        just = c("center", "bottom"), gp = gpar(cex = opt$`label-cex`, fontface = 3)
      )

      popViewport()
      popViewport()
      k <- k + 1
    }
  }
  popViewport()
  popViewport()
}

pt_list <- read_list(opt$`pt-list`)
mt_list <- read_list(opt$`mt-list`)
apt <- read_annot(opt$`annot-pt`, org = "pt")
amt <- read_annot(opt$`annot-mt`, org = "mt")

pdf(opt$out, width = opt$`page-width-in`, height = opt$`page-height-in`)
draw_page(pt_list, apt, opt$`pt-title`, org = "pt")
draw_page(mt_list, amt, opt$`mt-title`, org = "mt")
dev.off()
cat(opt$out, "\n")
