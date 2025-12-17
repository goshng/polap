#!/usr/bin/env Rscript
# scripts/polap-r-oatk-ptmt-sheet.R
# Version : v1.8.0  (2025-12-15)
# Author  : Sang Chul Choi (POLAP)
# License : GPL-3.0+
#
# v1.8.0 changes (relative to v1.7.1)
# - One page per species dataset
# v1.7.1 changes (relative to v1.7.0)
# - Keep ALL v1.7.0 pages (species × pattern) unchanged
# - ADD: one extra "summary" page per species that shows
#        tippo, himt, oatk-x-10, oatk-x-20, oatk-x-30 in a single page

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
    help = "CSV: species,code2,pattern,cov,mt_png,pt_png"
  ),
  make_option(
    c("--out"),
    type = "character",
    help = "Output PDF file"
  ),
  make_option(
    c("--title"),
    type    = "character",
    default = "OATK mitogenome assembly graphs (per species, per pattern)"
  ),
  make_option(c("--page-width-in"), type = "double", default = 11),
  make_option(c("--page-height-in"), type = "double", default = 17),
  make_option(c("--margin-in"), type = "double", default = 0.5),
  make_option(c("--gap-mm"), type = "double", default = 4.0),
  make_option(c("--label-cex"), type = "double", default = 0.9),
  make_option(c("--title-cex"), type = "double", default = 1.1),
  make_option(
    c("--rows-per-page"),
    type    = "integer",
    default = 5,
    help    = "Rows of tiles per page (default: 5)"
  ),
  make_option(
    c("--columns-per-page"),
    type    = "integer",
    default = 8,
    help    = "Columns of tiles per page (default: 8)"
  ),
  make_option(
    c("--start-page"),
    type    = "integer",
    default = 1,
    help    = "Starting page number for footer (default: 1)"
  ),
  # compatibility shim, ignored (v1.7.0 behavior)
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
  # Keep v1.7.1
  # per-pattern: current output (summary + pattern pages)
  # per-dataset: new output (one dataset page with 40 tiles)
  make_option(
    c("--mode"),
    type = "character",
    default = "per-dataset",
    help = "Page mode: per-pattern or per-dataset (default)"
  ),
  # MT/PT selection
  make_option(
    c("--mt-only"),
    action  = "store_true",
    default = FALSE,
    help    = "Show MT graphs only (mt_png column)"
  ),
  make_option(
    c("--pt-only"),
    action  = "store_true",
    default = FALSE,
    help    = "Show PT graphs only (pt_png column)"
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

# infer "summary tool" label from a single row
infer_summary_tool <- function(pattern, cov, mt_png, pt_png) {
  pat <- ifelse(is.na(pattern), "", as.character(pattern))
  cov_i <- suppressWarnings(as.integer(to_num(cov)))
  p_any <- ifelse(!is.na(mt_png) && nzchar(mt_png), mt_png, pt_png)

  if (!is.na(p_any) && nzchar(p_any)) {
    if (grepl("/himt/", p_any, fixed = TRUE) || grepl("himt_", basename(p_any), fixed = TRUE)) {
      return("himt")
    }
    if (grepl("/tippo/", p_any, fixed = TRUE)) {
      return("tippo")
    }
  }

  if (pat == "x" && !is.na(cov_i) && cov_i %in% c(10, 20, 30)) {
    return(sprintf("oatk-x-%d", cov_i))
  }

  NA_character_
}

# build a 5-tile "summary" slice (one row per target tool)
make_summary_slice <- function(df_sp_all) {
  targets <- c("oatk-x-10", "oatk-x-20", "oatk-x-30", "himt", "tippo")

  df_sp_all2 <- df_sp_all %>%
    mutate(.sum_tool = mapply(infer_summary_tool, pattern, cov, mt_png, pt_png, SIMPLIFY = TRUE))

  # pick a representative row for each target tool; if missing, create blank row
  out_rows <- lapply(targets, function(tt) {
    hit <- df_sp_all2 %>% filter(.sum_tool == tt)
    if (nrow(hit) > 0) {
      # choose the first (cov-sorted already in caller)
      r <- hit[1, , drop = FALSE]
      r$.tile_label <- tt
      return(r)
    } else {
      # create a blank placeholder tile
      r <- df_sp_all2[1, , drop = FALSE]
      r$pattern <- NA
      r$cov <- NA
      r$completeness_prop <- NA
      r$mt_png <- NA_character_
      r$pt_png <- NA_character_
      r$.tile_label <- tt
      return(r)
    }
  })

  bind_rows(out_rows)
}

# ------------------------------------------------------------------#
# Data: main CSV
# ------------------------------------------------------------------#
df <- suppressMessages(read_csv(
  opt$list,
  show_col_types = FALSE
))

# accept both column names (your CSV uses completeness_prop)
if (!("completeness_prop" %in% names(df)) && ("geneset_completeness_prop" %in% names(df))) {
  df$completeness_prop <- df$geneset_completeness_prop
}

required_cols <- c(
  "species",
  "code2",
  "pattern",
  "cov",
  "mt_png",
  "pt_png",
  "completeness_prop"
)
# Add tippo and himt.
if ("tool" %in% names(df)) required_cols <- c(required_cols, "tool")

missing_cols <- setdiff(required_cols, names(df))
if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in --list CSV: ",
    paste(missing_cols, collapse = ", ")
  )
}

df <- df[, required_cols] %>%
  mutate(
    cov = to_num(cov),
    completeness_prop = to_num(completeness_prop),
    pattern = factor(pattern)
  )

pattern_levels <- levels(df$pattern)
pattern_levels <- pattern_levels[!is.na(pattern_levels)]

if (length(pattern_levels) == 0) {
  stop("No patterns found in CSV.")
}

# the 40-case order once
# the target tile order that a dataset page must follow.
cov_all <- c(2, 5, 8, 11, 20, 30, 40, 50)
cov_xx <- c(11, 20, 30, 40, 50)

cases_40 <- c(
  sprintf("t-h-%02d", cov_all),
  sprintf("h-t-%02d", cov_all),
  sprintf("t-x-%02d", cov_all),
  sprintf("x-h-%02d", cov_all),
  sprintf("x-x-%02d", cov_xx),
  "tippo", "himt", "x-30"
)

# which PNG column to use
image_col <- if (opt$`mt-only`) {
  "mt_png"
} else if (opt$`pt-only`) {
  "pt_png"
} else {
  "mt_png"
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
# Build a case_id for every row (v1.8.0)
# Turns a row into a single canonical slot in the 40-case page.
# ------------------------------------------------------------------#
make_case_id <- function(tool, pattern, cov) {
  tool <- as.character(tool)
  pat <- ifelse(is.na(pattern), "", as.character(pattern))
  cv_i <- suppressWarnings(as.integer(round(as.numeric(cov))))

  if (tool == "tippo") {
    return("tippo")
  }
  if (tool == "himt") {
    return("himt")
  }

  # patterned oatk
  if (tool == "oatk" && pat %in% c("t-h", "h-t", "t-x", "x-h", "x-x") && !is.na(cv_i)) {
    return(sprintf("%s-%02d", pat, cv_i))
  }

  # baseline oatk x-30
  if (tool == "oatk" && pat %in% c("x", "") && !is.na(cv_i) && cv_i == 30L) {
    return("x-30")
  }

  NA_character_
}

# ------------------------------------------------------------------#
# New builder: one page per dataset with all 40 slots (v1.8.0)
# Ensures exactly 40 rows in the correct order
# Missing cases become NA -> your existing code prints No assembly
# Provides .tile_label so you can label tiles by t-h-02, etc.
# ------------------------------------------------------------------#
make_dataset_slice <- function(df_sp_all, cases_40) {
  # df_sp_all is already one species+code2
  d <- df_sp_all %>%
    mutate(
      tool = if ("tool" %in% names(.)) as.character(tool) else "oatk", # if your CSV lacks tool, set it earlier
      case_id = mapply(make_case_id, tool, pattern, cov, SIMPLIFY = TRUE)
    ) %>%
    filter(!is.na(case_id)) %>%
    group_by(case_id) %>%
    slice(1) %>% # if duplicates exist, keep first
    ungroup()

  # skeleton with all cases (keeps order and forces missing tiles)
  skeleton <- data.frame(case_id = cases_40, stringsAsFactors = FALSE)

  out <- skeleton %>%
    left_join(d, by = "case_id")
  sp_name0 <- df_sp_all$species[1]
  sp_code20 <- df_sp_all$code2[1]
  out$species <- sp_name0
  out$code2 <- sp_code20


  # Use case_id as the tile label for non-summary pages
  out$.tile_label <- out$case_id

  out
}



# ------------------------------------------------------------------#
# Page drawing (v1.7.0 + summary-mode extension)
# ------------------------------------------------------------------#
draw_pattern_page <- function(df_pat, pattern_label,
                              global_page_index, n_pages_global,
                              species_codes_df, subpage_index, n_subpages,
                              image_col,
                              .is_summary = FALSE) {
  n_tiles <- nrow(df_pat)

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
  # main_title <- opt$title
  # grid.text(
  #   main_title,
  #   x = 0.5,
  #   y = unit(1, "npc"),
  #   just = c("center", "top"),
  #   gp = gpar(fontface = "bold", cex = opt$`title-cex`)
  # )
  #
  # df_sp_uniq <- df_pat %>%
  #   distinct(species, code2)
  #
  # if (nrow(df_sp_uniq) == 1L) {
  #   sp_name <- df_sp_uniq$species[1]
  #   sp_code2 <- df_sp_uniq$code2[1]
  #   sp_pretty <- gsub("_", " ", sp_name, fixed = TRUE)
  #   species_hdr <- sprintf("%s (%s)", sp_code2, sp_pretty)
  # } else {
  #   species_hdr <- "Multiple species"
  # }
  #
  # sel_str <- if (image_col == "mt_png") "MT" else "PT"
  #
  # subtitle <- if (.is_summary) {
  #   sprintf("Species: %s  |  Summary: oatk-x-10/20/30 + himt + tippo (%s)", species_hdr, sel_str)
  # } else if (n_subpages > 1L) {
  #   sprintf(
  #     "Species: %s  |  Pattern: %s (%s, page %d/%d for this species+pattern)",
  #     species_hdr, pattern_label, sel_str, subpage_index, n_subpages
  #   )
  # } else {
  #   sprintf(
  #     "Species: %s  |  Pattern: %s (%s)",
  #     species_hdr, pattern_label, sel_str
  #   )
  # }
  #
  # grid.text(
  #   subtitle,
  #   x = 0.5,
  #   y = unit(1, "npc") - unit(1.4, "lines"),
  #   just = c("center", "top"),
  #   gp = gpar(cex = 0.9)
  # )

  # ----- Build dataset label (single species expected for all-cases) -----
  df_sp_uniq <- df_pat %>% distinct(species, code2)

  sp_name <- if (nrow(df_sp_uniq) >= 1L) df_sp_uniq$species[1] else ""
  sp_code2 <- if (nrow(df_sp_uniq) >= 1L) df_sp_uniq$code2[1] else ""
  sp_pretty <- gsub("_", " ", sp_name, fixed = TRUE)

  # Line 1: title + dataset
  line1 <- sprintf("%s - %s %s", opt$title, sp_code2, sp_pretty)

  # Line 2: fixed Methods hint
  line2 <- "Tiara for t, HiMT filter for h, no filtering for x, and the syncmer coverage threshold followed by gene-set completeness"

  left_x <- unit(0.01, "npc") # small inset from left margin

  grid.text(
    line1,
    x = left_x,
    y = unit(1, "npc"),
    just = c("left", "top"),
    gp = gpar(fontface = "bold", cex = opt$`title-cex`)
  )

  grid.text(
    line2,
    x = left_x,
    y = unit(1, "npc") - unit(1.2, "lines"),
    just = c("left", "top"),
    gp = gpar(cex = 0.85)
  )

  # Title: center align
  # grid.text(
  #   line1,
  #   x = 0.5,
  #   y = unit(1, "npc"),
  #   just = c("left", "top"),
  #   gp = gpar(fontface = "bold", cex = opt$`title-cex`)
  # )
  #
  # grid.text(
  #   line2,
  #   x = 0.5,
  #   y = unit(1, "npc") - unit(1.2, "lines"),
  #   just = c("left", "top"),
  #   gp = gpar(cex = 0.85)
  # )


  # ----------------------- Species legend (per page) --------------- #
  # legend_y <- unit(1, "npc") - unit(2.5, "lines")
  # # content_top <- unit(1, "npc") - unit(3.2, "lines")
  # # We use only 2 header lines now; start content higher
  # content_top <- unit(1, "npc") - unit(2.2, "lines")
  #
  #
  # if (!is.null(species_codes_df)) {
  #   df_uniq <- df_pat %>%
  #     distinct(species, code2) %>%
  #     filter(!is.na(code2), nzchar(code2), !is.na(species), nzchar(species))
  #
  #   legend_df <- df_uniq %>%
  #     left_join(species_codes_df, by = c("species" = "species"))
  #
  #   legend_df <- legend_df %>%
  #     mutate(
  #       code_label    = if_else(!is.na(code), code, code2),
  #       species_label = gsub("_", " ", species, fixed = TRUE),
  #       legend_label  = paste(code_label, species_label)
  #     )
  #
  #   legend_text <- paste(unique(legend_df$legend_label), collapse = "   ")
  #
  #   if (nzchar(legend_text)) {
  #     grid.text(
  #       legend_text,
  #       x = 0.5,
  #       y = legend_y,
  #       just = c("center", "top"),
  #       gp = gpar(cex = 0.85)
  #     )
  #   }
  # } else {
  #   df_uniq <- df_pat %>%
  #     distinct(species, code2) %>%
  #     filter(!is.na(code2), nzchar(code2), !is.na(species), nzchar(species)) %>%
  #     mutate(
  #       species_label = gsub("_", " ", species, fixed = TRUE),
  #       legend_label  = paste(code2, species_label)
  #     )
  #   legend_text <- paste(unique(df_uniq$legend_label), collapse = "   ")
  #   if (nzchar(legend_text)) {
  #     grid.text(
  #       legend_text,
  #       x = 0.5,
  #       y = legend_y,
  #       just = c("center", "top"),
  #       gp = gpar(cex = 0.85)
  #     )
  #   }
  # }

  # ----------------------- (Disabled) Species legend --------------- #
  # We keep only the 2-line header (title + Methods note).
  # So do NOT draw the per-page species legend line.
  content_top <- unit(1, "npc") - unit(2.2, "lines")


  # ----------------------- Body layout --------------------------- #
  footer_height <- unit(1.2, "lines")
  header_body_gap <- unit(1.0, "lines")
  body_height <- content_top - footer_height - header_body_gap

  body_vp <- viewport(
    y      = footer_height,
    height = body_height,
    width  = unit(1, "npc"),
    just   = c("center", "bottom")
  )
  pushViewport(body_vp)

  # For summary page: force a single-row 5-column grid
  # if (.is_summary) {
  #   sum_rows <- 1L
  #   sum_cols <- 5L
  #   pushViewport(viewport(
  #     layout = grid.layout(
  #       nrow    = sum_rows,
  #       ncol    = sum_cols,
  #       widths  = unit(rep(1, sum_cols), "null"),
  #       heights = unit(rep(1, sum_rows), "null")
  #     )
  #   ))
  #   max_tiles <- sum_rows * sum_cols
  #   local_rows <- sum_rows
  #   local_cols <- sum_cols
  if (.is_summary) {
    sum_cols <- 4L
    sum_rows <- ceiling(n_tiles / sum_cols)

    pushViewport(viewport(
      layout = grid.layout(
        nrow = sum_rows,
        ncol = sum_cols,
        widths = unit(rep(1, sum_cols), "null"),
        heights = unit(rep(1, sum_rows), "null")
      )
    ))

    max_tiles <- sum_rows * sum_cols
    local_rows <- sum_rows
    local_cols <- sum_cols

    # (keep your summary image sizing numbers or adjust slightly)
    img_height_frac <- 0.48
    img_y_center <- 0.42
  } else {
    pushViewport(viewport(
      layout = grid.layout(
        nrow    = rows_per_page,
        ncol    = cols_per_page,
        widths  = unit(rep(1, cols_per_page), "null"),
        heights = unit(rep(1, rows_per_page), "null")
      )
    ))
    max_tiles <- rows_per_page * cols_per_page
    local_rows <- rows_per_page
    local_cols <- cols_per_page

    # img_height_frac <- 0.75
    # img_y_center <- 0.45
  }

  for (k in seq_len(max_tiles)) {
    row_idx <- ((k - 1L) %/% local_cols) + 1L
    col_idx <- ((k - 1L) %% local_cols) + 1L
    if (row_idx > local_rows) next

    pushViewport(viewport(
      layout.pos.row = row_idx,
      layout.pos.col = col_idx
    ))

    pushViewport(viewport(
      x = 0.5, y = 0.5,
      width = unit(1, "npc") - 2 * gapx,
      height = unit(1, "npc") - 2 * gapy,
      # clip = "off"
      clip = "on"
    ))

    if (k <= n_tiles) {
      # Subtitle for each subfigure: p=t-h c=20
      # cov <- df_pat$cov[k]
      # pat_k <- as.character(df_pat$pattern[k])
      # img_path <- df_pat[[image_col]][k]
      #
      # # label per subfigure:
      # if (.is_summary) {
      #   row_label <- df_pat$.tile_label[k]
      # } else {
      #   row_label <- sprintf("p=%s c=%s", pat_k, cov)
      # }

      # Subtitle for each subfigure: t-h 20
      cov <- df_pat$cov[k]
      pat_k <- as.character(df_pat$pattern[k])
      img_path <- df_pat[[image_col]][k]

      # completeness (geneset) value if available
      comp_raw <- NA
      if ("completeness_prop" %in% names(df_pat)) comp_raw <- df_pat$completeness_prop[k]
      if (is.na(comp_raw) && "geneset_completeness_prop" %in% names(df_pat)) comp_raw <- df_pat$geneset_completeness_prop[k]

      comp_num <- suppressWarnings(as.numeric(comp_raw))
      comp_str <- if (!is.na(comp_num)) sprintf("%.2f", comp_num) else ""

      # v1.7.1
      # if (.is_summary) {
      #   # e.g. oatk-x-10-0.95  (or just oatk-x-10 if comp missing)
      #   base_lab <- df_pat$.tile_label[k]
      #   row_label <- if (nzchar(comp_str)) paste0(base_lab, "-", comp_str) else base_lab
      # } else {
      #   # e.g. t-h-20-0.95 (two decimals)
      #   cov_i <- suppressWarnings(as.integer(round(as.numeric(cov))))
      #   cov_str <- if (!is.na(cov_i)) as.character(cov_i) else as.character(cov)
      #   base_lab <- paste0(pat_k, "-", cov_str)
      #   row_label <- if (nzchar(comp_str)) paste0(base_lab, "-", comp_str) else base_lab
      # }

      # v1.8.0
      # Summary pages already set .tile_label (oatk-x-10, himt, tippo) -> used automatically.
      # Your future per-dataset pages will also set .tile_label (t-h-02, ..., x-30) -> used automatically.
      # Old per-pattern pages don’t have .tile_label -> it falls back to pattern-cov formatting.
      # Completeness suffix (-0.95) still appends.
      # v1.8.0
      # If a tile label was precomputed (summary page OR per-dataset page), use it.
      # Otherwise fall back to pattern+cov.
      if (".tile_label" %in% names(df_pat) && !is.na(df_pat$.tile_label[k]) && nzchar(as.character(df_pat$.tile_label[k]))) {
        base_lab <- as.character(df_pat$.tile_label[k])
      } else {
        cov_i <- suppressWarnings(as.integer(round(as.numeric(cov))))
        cov_str <- if (!is.na(cov_i)) sprintf("%02d", cov_i) else as.character(cov) # keep 02-style
        base_lab <- paste0(pat_k, "-", cov_str)
      }

      row_label <- if (nzchar(comp_str)) paste0(base_lab, "-", comp_str) else base_lab


      # Change scheme-cov font size
      label_cex_use <- opt$`label-cex`
      if (".tile_label" %in% names(df_pat)) {
        label_cex_use <- opt$`label-cex` * 0.85
      }
      grid.text(
        row_label,
        x = 0,
        y = 0.995,
        just = c("left", "top"),
        gp = gpar(cex = label_cex_use, fontface = 3)
      )

      img_height_frac <- 0.75
      img_y_center <- 0.45

      if (!is.na(img_path) && nzchar(img_path) && file.exists(img_path)) {
        im <- try(readPNG(img_path), silent = TRUE)

        if (!inherits(im, "try-error")) {
          if (.is_summary) {
            # Preserve aspect ratio by leaving height = NULL
            # (grid docs: if both width & height are set, image is likely distorted)
            grid.raster(
              im,
              width = unit(0.95, "npc"),
              height = NULL,
              y = unit(0.45, "npc"),
              interpolate = TRUE
            )
          } else {
            # keep your original behavior for non-summary pages
            grid.raster(
              im,
              width = unit(0.95, "npc"),
              height = unit(img_height_frac, "npc"),
              y = unit(img_y_center, "npc"),
              interpolate = TRUE
            )
          }
        } else {
          grid.text("x",
            gp = gpar(col = "red3", cex = 2.0, fontface = "bold"),
            y  = unit(img_y_center, "npc")
          )
        }
      } else {
        grid.text("x",
          gp = gpar(col = "red3", cex = 2.0, fontface = "bold"),
          y  = unit(img_y_center, "npc")
        )
      }
    }

    popViewport()
    popViewport()
  }

  popViewport() # layout
  popViewport() # body_vp

  # ----------------------- Footer (page number) ------------------- #
  global_page_no <- opt$`start-page` + global_page_index - 1L
  last_global_pg <- opt$`start-page` + n_pages_global - 1L
  footer_l <- sprintf("page %d/%d", global_page_no, last_global_pg)

  grid.text(
    footer_l,
    x = 0.5,
    y = unit(0.3, "lines"),
    just = c("center", "bottom"),
    gp = gpar(cex = 0.8)
  )

  popViewport()
}

# ------------------------------------------------------------------#
# Build page list (v1.8.0 and v1.7.1)
# If --mode per-pattern, you keep your working v1.7.1 output.
# Otherwise, you get one page per dataset with 40 tiles.
# ------------------------------------------------------------------#
page_specs <- list()

species_list <- df %>%
  distinct(species, code2) %>%
  arrange(species, code2)

for (sp_i in seq_len(nrow(species_list))) {
  sp_name <- species_list$species[sp_i]
  sp_code2 <- species_list$code2[sp_i]

  df_sp_all <- df %>%
    filter(species == sp_name, code2 == sp_code2) %>%
    arrange(cov)

  if (nrow(df_sp_all) == 0L) next

  if (opt$mode == "per-dataset") {
    df_page <- make_dataset_slice(df_sp_all, cases_40)

    page_specs[[length(page_specs) + 1L]] <- list(
      species     = sp_name,
      code2       = sp_code2,
      pattern     = "all-cases",
      df          = df_page,
      subpage_idx = 1L,
      n_subpages  = 1L,
      is_summary  = FALSE
    )
  } else {
    # ---- current v1.7.1 behavior ----
    df_sum <- make_summary_slice(df_sp_all)
    page_specs[[length(page_specs) + 1L]] <- list(
      species     = sp_name,
      code2       = sp_code2,
      pattern     = "summary",
      df          = df_sum,
      subpage_idx = 1L,
      n_subpages  = 1L,
      is_summary  = TRUE
    )

    for (pat in pattern_levels) {
      df_sp_pat_all <- df_sp_all %>%
        filter(pattern == pat) %>%
        arrange(cov)

      if (nrow(df_sp_pat_all) == 0L) next

      n_tiles_pat <- nrow(df_sp_pat_all)
      n_chunks <- ceiling(n_tiles_pat / tiles_per_page)

      for (chunk_i in seq_len(n_chunks)) {
        start_i <- (chunk_i - 1L) * tiles_per_page + 1L
        end_i <- min(start_i + tiles_per_page - 1L, n_tiles_pat)
        df_slice <- df_sp_pat_all[start_i:end_i, , drop = FALSE]

        page_specs[[length(page_specs) + 1L]] <- list(
          species     = sp_name,
          code2       = sp_code2,
          pattern     = as.character(pat),
          df          = df_slice,
          subpage_idx = chunk_i,
          n_subpages  = n_chunks,
          is_summary  = FALSE
        )
      }
    }
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
    df_pat             = spec$df,
    pattern_label      = spec$pattern,
    global_page_index  = i,
    n_pages_global     = n_pages_global,
    species_codes_df   = species_codes_df,
    subpage_index      = spec$subpage_idx,
    n_subpages         = spec$n_subpages,
    image_col          = image_col,
    .is_summary        = isTRUE(spec$is_summary)
  )
}

dev.off()
cat(opt$out, "\n")
