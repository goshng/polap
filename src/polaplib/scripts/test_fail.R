#!/usr/bin/env Rscript
# ##############################################################################
# This file is part of polap. GPLv3+ (see top-level license)
# ##############################################################################
# test_fail.R
# Intentional R failure with backtrace printing and exit code control.
# Usage: test_fail.R [--message STR] [--exit-code N] [INFILE [OUTFILE]]

args <- commandArgs(trailingOnly = TRUE)

message <- "intentional test failure (R)"
exit_code <- 3L
infile <- ""
outfile <- ""

# Minimal arg parsing
i <- 1L
while (i <= length(args)) {
  a <- args[[i]]
  if (a == "--message" && i + 1L <= length(args)) {
    message <- args[[i + 1L]]
    i <- i + 2L
  } else if (a == "--exit-code" && i + 1L <= length(args)) {
    exit_code <- as.integer(args[[i + 1L]])
    i <- i + 2L
  } else if (a == "-h" || a == "--help") {
    cat("Usage: test_fail.R [--message STR] [--exit-code N] [INFILE [OUTFILE]]\n")
    quit(save = "no", status = 0L)
  } else if (substr(a, 1L, 2L) == "--") {
    stop(sprintf("Unknown option: %s", a))
  } else {
    # positional
    if (identical(infile, "")) infile <- a else if (identical(outfile, "")) outfile <- a
    i <- i + 1L
  }
}

# Optional write
if (nzchar(outfile)) {
  con <- file(outfile, open = "w")
  writeLines("test_fail.R writing before failure", con)
  close(con)
}

# Ensure we get a printed backtrace and controlled exit status in non-interactive runs
mk_handler <- function(code) {
  function() {
    # Base traceback (last error). Depth 2 prints calls; adjust as desired.
    utils::traceback(2)
    # Exit with the requested status
    quit(save = "no", status = code)
  }
}
options(error = mk_handler(exit_code))

third <- function() {
  stop(message, call. = TRUE)
}
second <- function() third()
first <- function() second()

# Emit a little context and then fail
message(sprintf("R test_fail: message='%s' code=%d infile='%s' outfile='%s'", message, exit_code, infile, outfile))
first()
