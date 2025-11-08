#!/usr/bin/env Rscript
# polaplib/polap-r-test-crash.R
# Version: v0.1.1  (optparse; stack ON by default; quiet; robust)

suppressPackageStartupMessages(library(optparse))
options(show.error.messages = FALSE) # our handler controls error output

# --- CLI (optparse) -----------------------------------------------------------
# NOTE:
#  - We use *one* option for exit code: -c / --exit-code (no extra alias) to avoid
#    optparse internals confusion.
#  - We set dest names without hyphens (e.g., no_stack) so we can access as opts$no_stack.
option_list <- list(
  make_option(c("-m", "--message"),
    type    = "character",
    default = "intentional test failure (R)",
    help    = "failure message [default %default]"
  ),
  make_option(c("-c", "--exit-code"),
    type    = "integer",
    default = 7,
    dest    = "exit_code",
    help    = "non-zero exit code [default %default]"
  ),
  make_option(c("--in"),
    type    = "character",
    default = "",
    dest    = "infile",
    help    = "optional input file"
  ),
  make_option(c("--out"),
    type    = "character",
    default = "",
    dest    = "outfile",
    help    = "optional output file (written pre-crash)"
  ),
  make_option(c("--no-stack"),
    action  = "store_true",
    default = FALSE,
    dest    = "no_stack",
    help    = "suppress short stack (STACK IS ON BY DEFAULT)"
  )
)

parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
opts <- parse_args(parser) # standard: returns a named list of options

# Normalize/guard options
msg <- if (!is.null(opts$message)) opts$message else ""
ecode <- if (!is.null(opts$exit_code)) opts$exit_code else 7L
if (is.na(ecode) || ecode == 0L) ecode <- 7L
infile <- if (!is.null(opts$infile)) opts$infile else ""
outfile <- if (!is.null(opts$outfile)) opts$outfile else ""
# STACK ON by default; user can turn it off with --no-stack
wantStack <- !isTRUE(opts$no_stack)

# --- Error handlers (quiet by default; stack prints short call chain) ---------
quiet_handler <- function() {
  quit(save = "no", status = as.integer(ecode))
}

stack_handler <- local({
  print_stack <- function() {
    calls <- sys.calls()
    nc <- length(calls)
    if (!nc) {
      cat("Stack (empty)\n")
      return()
    }
    cat("Stack (oldest → newest):\n")
    for (i in seq_len(nc)) {
      cl <- calls[[i]]
      sr <- attr(cl, "srcref")
      txt <- paste(deparse(cl), collapse = " ")
      if (!is.null(sr)) {
        sf <- attr(sr, "srcfile")
        fn <- if (!is.null(sf) && !is.null(sf$filename)) basename(sf$filename) else "<src>"
        cat(sprintf(
          "#%d %s:%d:%d  %s\n",
          i, fn, as.integer(sr[[1]]), as.integer(sr[[2]]), txt
        ))
      } else {
        cat(sprintf("#%d %s\n", i, txt))
      }
    }
  }
  function() {
    try(print_stack(), silent = TRUE)
    quit(save = "no", status = as.integer(ecode))
  }
})

options(error = if (wantStack) stack_handler else quiet_handler)

# --- Context & optional pre-crash side effect ---------------------------------
cat(sprintf(
  "R test_crash: message='%s' code=%d infile='%s' outfile='%s'\n",
  msg, as.integer(ecode), infile, outfile
))

if (nzchar(outfile)) {
  con <- file(outfile, open = "w")
  writeLines("polap-r-test-crash.R wrote before failing", con)
  close(con)
}

# --- Nested failure to exercise the handler -----------------------------------
third <- function(m, c) stop(m, call. = TRUE)
second <- function(m, c) third(m, c)
first <- function(m, c) second(m, c)

first(msg, as.integer(ecode))
