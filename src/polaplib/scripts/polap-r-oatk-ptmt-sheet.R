#!/usr/bin/env Rscript
# scripts/polap-r-oatk-ptmt-sheet.R
# Version : v1.6.0  (2025-12-08)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# Draw OATK mitogenome (or plastid) assembly graphs per pattern, laid out
# on a rows-per-page × columns-per-page grid.
#
# Each CSV row corresponds to ONE tile (subfigure):
#   species, code2, pattern, cov, mito_png, plastid_png
#
# You choose which PNG to show:
#   --mt-only   → use mito_png
#   --pt-only   → use plastid_png
#   (default: mito_png)
#
# Grid:
#   --rows-per-page R
#   --columns-per-page C
#   → R × C tiles per page
#
# Pages are grouped by pattern; if a pattern has more than R*C rows, it
# gets multiple pages: “Pattern: h-t (page 2/3 for this pattern)” etc.
#
# Species-codes file (--species-codes):
#   space-delimited, columns:
#       code species
#   e.g.,
#       Au01 Aegilops_umbellulata
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
  make_option(
    c("--list"),
    type = "character",
    help = "CSV: species,code2,pattern,cov,mito_png,plastid_png"
  ),
  make_option(
    c("--out"),
    type = "character",
    help = "Output PDF file"
  ),
  make_option(
    c("--title"),
    type    = "character",
    default = "OATK mitogenome assembly graphs (per pattern)"
  ),
  make_option(c("--page-width-in"),  type = "double", default = 11),
  make_option(c("--page-height-in"), type = "double", default = 17),
  make_option(c("--margin-in"),      type = "double", default = 0.5),
  make_option(c("--gap-mm"),         type = "double", default = 4.0),
  make_option(c("--label-cex"),      type = "double", default = 0.9),
  make_option(c("--title-cex"),      type = "double", default = 1.1),
  make_option(
    c("--rows-per-page"),
    type    = "integer",
    default = 8,
    help    = "Rows of tiles per page (default: 8)"
  ),
  make_option(
    c("--columns-per-page"),
    type    = "integer",
    default = 6,
    help    = "Columns of tiles per page (default: 6)"
  ),
  make_option(
    c("--start-page"),
    type    = "integer",
    default = 1,
    help    = "Starting page number for footer (default: 1)"
  ),
  # compatibility shim, ignored
  make_option(
    c("--base-dir"),
    type    = "character",
    default = ".",
    help    = "(ignored) kept for backward compatibility with old Make rules"
  ),
  # species-codes mapping file
  make_option(
    c("--species-codes"),
    type    = "character",
    default = NULL,
    help    = "Species-codes mapping file (space-delimited: code species)"
  ),
  # MT/PT selection
  make_option(
    c("--mt-only"),
    action  = "store_true",
    default = FALSE,
    help    = "Show MT graphs only (mito_png column)"
  ),
  make_option(
    c("--pt-only"),
    action  = "store_true",
    default = FALSE,
    help    = "Show PT graphs only (plastid_png column)"
  )
)) |>
  parse_args()

if (is.null(opt$list) || is.null(opt$out)) {
  stop("Need --list and --out")
}

if (opt$`mt-only` && opt$`pt-only`) {
  stop("Cannot use --mt-only and --pt-only together.")
}

# ------------------------------------------------------------------#
# Helpers
# ------------------------------------------------------------------#
to_num <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }
  suppressWarnings(readr::parse_number(as.character(x)))
}

# ------------------------------------------------------------------#
# Data: main CSV
# ------------------------------------------------------------------#
df <- suppressMessages(read_csv(
  opt$list,
  show_col_types = FALSE
))

required_cols <- c("species", "code2", "pattern", "cov", "mito_png", "plastid_png")
missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns in --list CSV: ",
       paste(missing_cols, collapse = ", "))
}

df <- df[, required_cols] %>%
  mutate(
    cov     = to_num(cov),
    pattern = factor(pattern)
  )

patterns <- levels(df$pattern)
patterns <- patterns[!is.na(patterns)]

if (length(patterns) == 0) {
  stop("No patterns found in CSV.")
}

# which PNG column to use
image_col <- if (opt$`mt-only`) {
  "mito_png"
} else if (opt$`pt-only`) {
  "plastid_png"
} else {
  # default: MT (this is a mitogenome sheet)
  "mito_png"
}

gapx <- unit(opt$`gap-mm`, "mm")
gapy <- unit(opt$`gap-mm`, "mm")

rows_per_page <- max(1L, as.integer(opt$`rows-per-page`))
cols_per_page <- max(1L, as.integer(opt$`columns-per-page`))
tiles_per_page <- rows_per_page * cols_per_page

# ------------------------------------------------------------------#
# Data: species-codes mapping (optional)
# ------------------------------------------------------------------#
species_codes_df <- NULL
if (!is.null(opt$`species-codes`) && nzchar(opt$`species-codes`) &&
    file.exists(opt$`species-codes`)) {
  species_codes_df <- suppressMessages(
    read_table(
      opt$`species-codes`,
      col_names = c("code", "species"),
      comment = "#",
      col_types = cols(
        code    = col_character(),
        species = col_character()
      )
    )
  )
}

# ------------------------------------------------------------------#
# Page drawing
# ------------------------------------------------------------------#
draw_pattern_page <- function(df_pat, pattern_label,
                              global_page_index, n_pages_global,
                              species_codes_df, subpage_index, n_subpages,
                              image_col) {
  n_tiles <- nrow(df_pat)

  grid.newpage()

  # Outer viewport with margins
  m <- unit(opt$`margin-in`, "in")
  outer_vp <- viewport(
    x = 0.5, y = 0.5,
    width  = unit(1, "npc") - 2 * m,
    height = unit(1, "npc") - 2 * m
  )
  pushViewport(outer_vp)

  # ----------------------- Title & subtitle ----------------------- #
  main_title <- opt$title
  grid.text(
    main_title,
    x = 0.5,
    y = unit(1, "npc"),
    just = c("center", "top"),
    gp = gpar(fontface = "bold", cex = opt$`title-cex`)
  )

  sel_str <- if (image_col == "mito_png") {
    "MT"
  } else if (image_col == "plastid_png") {
    "PT"
  } else {
    "MT"
  }

  if (n_subpages > 1L) {
    subtitle <- sprintf("Pattern: %s (%s, page %d/%d for this pattern)",
                        pattern_label, sel_str, subpage_index, n_subpages)
  } else {
    subtitle <- sprintf("Pattern: %s (%s)", pattern_label, sel_str)
  }

  grid.text(
    subtitle,
    x = 0.5,
    y = unit(1, "npc") - unit(1.4, "lines"),
    just = c("center", "top"),
    gp = gpar(cex = 0.9)
  )

  # ----------------------- Species legend (per page) --------------- #
  legend_y    <- unit(1, "npc") - unit(2.5, "lines")
  content_top <- unit(1, "npc") - unit(3.2, "lines")

  if (!is.null(species_codes_df)) {
    df_uniq <- df_pat %>%
      distinct(species, code2)

    legend_df <- df_uniq %>%
      left_join(species_codes_df, by = c("species" = "species"))

    legend_df <- legend_df %>%
      mutate(
        code_label    = if_else(!is.na(code), code, code2),
        species_label = gsub("_", " ", species, fixed = TRUE),
        legend_label  = paste(code_label, species_label)
      )

    legend_text <- paste(unique(legend_df$legend_label), collapse = "   ")

    if (nzchar(legend_text)) {
      grid.text(
        legend_text,
        x = 0.5,
        y = legend_y,
        just = c("center", "top"),
        gp = gpar(cex = 0.85)
      )
    }
  } else {
    df_uniq <- df_pat %>%
      distinct(species, code2) %>%
      mutate(
        species_label = gsub("_", " ", species, fixed = TRUE),
        legend_label  = paste(code2, species_label)
      )
    legend_text <- paste(unique(df_uniq$legend_label), collapse = "   ")
    if (nzchar(legend_text)) {
      grid.text(
        legend_text,
        x = 0.5,
        y = legend_y,
        just = c("center", "top"),
        gp = gpar(cex = 0.85)
      )
    }
  }

  # ----------------------- Body layout --------------------------- #
  footer_height   <- unit(1.2, "lines")
  header_body_gap <- unit(1.0, "lines")
  body_height     <- content_top - footer_height - header_body_gap

  body_vp <- viewport(
    y      = footer_height,
    height = body_height,
    width  = unit(1, "npc"),
    just   = c("center", "bottom")
  )
  pushViewport(body_vp)

  pushViewport(viewport(
    layout = grid.layout(
      nrow    = rows_per_page,
      ncol    = cols_per_page,
      widths  = unit(rep(1, cols_per_page), "null"),
      heights = unit(rep(1, rows_per_page), "null")
    )
  ))

  # Place each tile in row-major order
  max_tiles <- rows_per_page * cols_per_page
  for (k in seq_len(max_tiles)) {
    row_idx <- ((k - 1L) %/% cols_per_page) + 1L
    col_idx <- ((k - 1L) %%  cols_per_page) + 1L

    if (row_idx > rows_per_page) {
      next
    }

    pushViewport(viewport(
      layout.pos.row = row_idx,
      layout.pos.col = col_idx
    ))

    pushViewport(viewport(
      x = 0.5, y = 0.5,
      width  = unit(1, "npc") - 2 * gapx,
      height = unit(1, "npc") - 2 * gapy,
      clip   = "off"
    ))

    if (k <= n_tiles) {
      sp  <- df_pat$species[k]
      cd2 <- df_pat$code2[k]
      cov <- df_pat$cov[k]
      img_path <- df_pat[[image_col]][k]

      # label per subfigure: coverage only, e.g. "c=1"
      row_label <- sprintf("c=%s", cov)

      grid.text(
        row_label,
        x = 0,
        y = 0.98,
        just = c("left", "top"),
        gp = gpar(cex = opt$`label-cex`, fontface = 3)
      )

      img_height_frac <- 0.75
      img_y_center    <- 0.45

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
    } else {
      # empty tile (optional: draw nothing)
    }

    popViewport() # inner
    popViewport() # cell
  }

  popViewport() # layout
  popViewport() # body_vp

  # ----------------------- Footer (page number) ------------------- #
  global_page_no  <- opt$`start-page` + global_page_index - 1L
  last_global_pg  <- opt$`start-page` + n_pages_global - 1L
  footer_l <- sprintf("page %d/%d", global_page_no, last_global_pg)

  grid.text(
    footer_l,
    x   = 0.5,
    y   = unit(0.3, "lines"),
    just = c("center", "bottom"),
    gp  = gpar(cex = 0.8)
  )

  popViewport() # outer_vp
}

# ------------------------------------------------------------------#
# Build page list (patterns × chunks of tiles_per_page)
# ------------------------------------------------------------------#
page_specs <- list()
for (pat in patterns) {
  df_pat_all <- df %>%
    filter(pattern == pat) %>%
    arrange(cov, species, code2)

  if (nrow(df_pat_all) == 0L) {
    next
  }

  n_tiles_pat <- nrow(df_pat_all)
  n_chunks    <- ceiling(n_tiles_pat / tiles_per_page)

  for (chunk_i in seq_len(n_chunks)) {
    start_i <- (chunk_i - 1L) * tiles_per_page + 1L
    end_i   <- min(start_i + tiles_per_page - 1L, n_tiles_pat)
    df_slice <- df_pat_all[start_i:end_i, , drop = FALSE]

    page_specs[[length(page_specs) + 1L]] <- list(
      pattern      = as.character(pat),
      df           = df_slice,
      subpage_idx  = chunk_i,
      n_subpages   = n_chunks
    )
  }
}

n_pages_global <- length(page_specs)
if (n_pages_global == 0L) {
  stop("No pages to draw (no data?).")
}

# ------------------------------------------------------------------#
# Generate PDF
# ------------------------------------------------------------------#
pdf(
  opt$out,
  width  = opt$`page-width-in`,
  height = opt$`page-height-in`
)

for (i in seq_along(page_specs)) {
  spec <- page_specs[[i]]
  draw_pattern_page(
    df_pat           = spec$df,
    pattern_label    = spec$pattern,
    global_page_index = i,
    n_pages_global   = n_pages_global,
    species_codes_df = species_codes_df,
    subpage_index    = spec$subpage_idx,
    n_subpages       = spec$n_subpages,
    image_col        = image_col
  )
}

dev.off()
cat(opt$out, "\n")
