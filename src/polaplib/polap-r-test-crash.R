#!/usr/bin/env Rscript
# polap-r-test-crash.R
# Version: v0.1.5
# Self-contained crash tester that prints file:line by re-reading its own file,
# so it also works under `Rscript --vanilla` where srcref isn't attached.

# --- find this script's path ---------------------------------------------------
.this_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  farg <- grep("^--file=", args, value = TRUE)
  if (length(farg)) {
    return(normalizePath(sub("^--file=", "", farg[1]), mustWork = FALSE))
  }
  # Fallbacks (sourced or unknown); last resort: basename in CWD
  if (!is.null(sys.frames()) && length(sys.frames())) {
    src <- attr(sys.frames()[[1]], "ofile")
    if (!is.null(src)) {
      return(normalizePath(src, mustWork = FALSE))
    }
  }
  # Last guess: assume we're running the file from its directory
  return(normalizePath("polap-r-test-crash.R", mustWork = FALSE))
}

# --- build a map: function name -> definition line in this file ---------------
.build_funline_map <- function(path) {
  ln <- tryCatch(readLines(path, warn = FALSE), error = function(...) character())
  if (!length(ln)) {
    return(list(lines = ln, map = list()))
  }
  # Simple regex: NAME <- function(
  rx <- "^\\s*([A-Za-z.][A-Za-z0-9._]*)\\s*<-\\s*function\\s*\\("
  m <- regexec(rx, ln, perl = TRUE)
  hits <- regmatches(ln, m)
  mp <- list()
  for (i in seq_along(hits)) {
    h <- hits[[i]]
    if (length(h) >= 2L) {
      name <- h[[2L]]
      if (!nzchar(name)) next
      # first definition wins
      if (is.null(mp[[name]])) mp[[name]] <- i
    }
  }
  list(lines = ln, map = mp)
}

# --- pretty backtrace using map (best effort) ----------------------------------
.print_bt <- function(e, funmap, script_path) {
  cat(sprintf("R-FAILED: %s\n", conditionMessage(e)), file = stderr())
  calls <- sys.calls()
  if (!length(calls)) {
    utils::traceback()
    return()
  }

  mp <- funmap$map
  depth <- 0L
  for (i in rev(seq_along(calls))) {
    cl <- calls[[i]]
    # function name
    fname <- tryCatch(
      {
        f <- cl[[1L]]
        if (is.symbol(f)) as.character(f) else as.character(f)[1]
      },
      error = function(...) "<unknown>"
    )

    # one-line call text
    call_txt <- tryCatch(paste(deparse(cl, width.cutoff = 200L), collapse = " "),
      error = function(...) "<call unavailable>"
    )

    if (!is.null(mp[[fname]])) {
      cat(sprintf("  #%d %s:%d: %s\n", depth, script_path, mp[[fname]], call_txt),
        file = stderr()
      )
    } else {
      cat(sprintf("  #%d %s\n", depth, call_txt), file = stderr())
    }
    depth <- depth + 1L
  }
}

# --- args & demo chain ---------------------------------------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  msg <- "hello from polap"
  code <- 1L
  i <- 1L
  while (i <= length(args)) {
    if (i < length(args) && args[i] == "--msg") {
      msg <- args[i + 1L]
      i <- i + 2L
      next
    }
    if (i < length(args) && args[i] == "--code") {
      code <- as.integer(args[i + 1L])
      i <- i + 2L
      next
    }
    i <- i + 1L
  }
  list(msg = msg, code = code)
}

leaf <- function(msg, code) stop(sprintf("demo crash (R): %s [code=%d]", msg, as.integer(code)))
mid <- function(msg, code) leaf(msg, code)
top <- function(msg, code) mid(msg, code)

main <- function() {
  a <- parse_args()
  top(a$msg, a$code)
}

# --- run with in-place handler and rethrow to exit non-zero --------------------
script_path <- .this_file()
funmap <- .build_funline_map(script_path)

status <- 0L
tryCatch(
  withCallingHandlers(
    main(),
    error = function(e) {
      .print_bt(e, funmap, script_path)
      # rethrow to propagate non-zero status
      stop(e)
    }
  ),
  error = function(e) {
    status <<- 10L
  }
)

quit(status = status, save = "no")
