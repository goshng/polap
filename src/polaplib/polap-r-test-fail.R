#!/usr/bin/env Rscript
# Version: v0.4.0
# Purpose: minimal crash-only R test using optparse + line-numbered stack report

suppressPackageStartupMessages({
  library(optparse)
})

## ── Minimal, clean stack: cut after the signaling frame ─────────────────────
options(keep.source = TRUE)

polap_error_handler <- local({
  # best effort: locate a file:line for a call
  frame_loc <- function(cl) {
    sr <- attr(cl, "srcref")
    if (!is.null(sr)) {
      srcfile  <- attr(sr, "srcfile")
      filename <- if (!is.null(srcfile) && !is.null(srcfile$filename))
                    basename(srcfile$filename) else ""
      start_ln <- as.integer(sr[[1]]); start_col <- as.integer(sr[[2]])
      return(list(file = filename, line = start_ln, col = start_col))
    }
    list(file = "", line = NA_integer_, col = NA_integer_)
  }

  # decide where to stop printing (at the error signal call)
  find_signal_cutoff <- function(calls) {
    # names that actually *signal* the error
    sig_syms <- c("stop", "stopifnot", "abort", "signalCondition",
                  "rlang::abort", "base::signalCondition")
    is_signal <- function(cl) {
      if (!is.call(cl)) return(FALSE)
      head <- cl[[1]]
      # as symbols / names / qualified calls
      nm <- tryCatch(as.character(head), error = function(e) character())
      any(nm %in% sig_syms)
    }
    idx <- which(vapply(calls, is_signal, logical(1)))
    if (length(idx)) {
      return(tail(idx, 1))  # last signal in the stack (closest to error point)
    }
    length(calls)  # fallback: print everything
  }

  print_stack <- function() {
    calls <- sys.calls()
    nc <- length(calls)
    if (!nc) { cat("No call stack available.\n"); return() }

    # cut stack AFTER the signaling frame so handler/tryCatch frames are gone
    cut <- find_signal_cutoff(calls)
    calls <- calls[seq_len(cut)]

    cat("Stack (oldest → newest):\n")
    for (i in seq_along(calls)) {
      cl <- calls[[i]]
      loc <- frame_loc(cl)
      call_txt <- paste(deparse(cl), collapse = " ")
      if (nzchar(loc$file) && !is.na(loc$line)) {
        cat(sprintf("#%d %s:%d:%d  %s\n",
                    i, loc$file, loc$line, ifelse(is.na(loc$col), 1L, loc$col), call_txt))
      } else {
        cat(sprintf("#%d %s\n", i, call_txt))
      }
    }
  }

  function(exit_code) {
    force(exit_code)
    function() {
      # Don’t introduce noisy frames; keep handler minimal & non-throwing
      tryCatch(print_stack(), error = function(e) {})
      quit(save = "no", status = exit_code)
    }
  }
})

## 3. Actual test logic -------------------------------------------------------
third  <- function(message_txt, exit_code) stop(message_txt, call. = TRUE)
second <- function(message_txt, exit_code) third(message_txt, exit_code)
run <- function(message_txt, exit_code, infile, outfile) {
  if (nzchar(outfile)) {
    con <- file(outfile, open = "w")
    writeLines("polap-r-test-fail.R wrote before failing", con)
    close(con)
  }
  cat(sprintf("R test_fail: message='%s' code=%d infile='%s' outfile='%s'\n",
              message_txt, exit_code, infile, outfile))
  first <- function() second(message_txt, exit_code)
  first()
}

## 4. Command handler (main) --------------------------------------------------
main <- function() {
  option_list <- list(
    make_option(c("-m", "--message"),   type = "character", default = "intentional test failure (R)",
                help = "failure message [default %default]"),
    make_option(c("-c", "--exit-code"), type = "integer", default = 7,
                help = "exit code [default %default]")
  )
  parser <- OptionParser(usage = "%prog [options] [INFILE [OUTFILE]]",
                         option_list = option_list,
                         description = "polap R test-fail script")

  argv <- parse_args(parser, positional_arguments = TRUE)
  infile  <- if (length(argv$args) >= 1) argv$args[[1]] else ""
  outfile <- if (length(argv$args) >= 2) argv$args[[2]] else ""

  # Install safe crash handler
  options(error = polap_error_handler(argv$options$`exit-code`))

  # Dispatch
  run(argv$options$message, argv$options$`exit-code`, infile, outfile)
}

## 5. Run only when executed directly ----------------------------------------
if (identical(environmentName(topenv()), "R_GlobalEnv")) {
  main()
}
