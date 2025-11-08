#!/usr/bin/env Rscript
# polaplib/scripts/test_crash.R
# Version: v0.2.4 (optparse; stack ON by default; quiet; robust to missing args)

suppressPackageStartupMessages(library(optparse))
options(keep.source = FALSE) # Rscript typically turns this off; OK for our use
options(show.error.messages = FALSE) # let our handler be the only error output

# ---- CLI: use optparse, stack ON by default, allow '-' in values -------------
option_list <- list(
  make_option(c("-m", "--message"),
    type = "character", default = "intentional test failure (R)",
    help = "failure message [default %default]"
  ),
  make_option(c("-c", "--exit-code"),
    type = "integer", default = 7, dest = "exit_code",
    help = "exit code [default %default]"
  ),
  make_option(c("--in"),
    type = "character", default = "", dest = "infile",
    help = "optional input file"
  ),
  make_option(c("--out"),
    type = "character", default = "", dest = "outfile",
    help = "optional output file (written before crash)"
  ),
  make_option(c("--no-stack"),
    action = "store_true", default = FALSE, dest = "no_stack",
    help = "suppress short stack (stack is ON by default)"
  )
)
parser <- OptionParser(usage = "%prog [options]", option_list = option_list)
argv <- parse_args(parser, positional_arguments = FALSE)

msg <- if (!is.null(argv$options$message)) argv$options$message else ""
ecode <- if (!is.null(argv$options$exit_code)) argv$options$exit_code else 7L
infile <- if (!is.null(argv$options$infile)) argv$options$infile else ""
outfile <- if (!is.null(argv$options$outfile)) argv$options$outfile else ""
wantStack <- !isTRUE(tryCatch(argv$options$`no_stack`, error = function(e) FALSE))

# ---- Error handler: quiet header only; optional short stack -------------------
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

# ---- Context + optional pre-crash side effect --------------------------------
cat(sprintf(
  "R test_crash: message='%s' code=%d infile='%s' outfile='%s'\n",
  msg, as.integer(ecode), infile, outfile
))

if (nzchar(outfile)) {
  con <- file(outfile, open = "w")
  writeLines("test_crash.R wrote before failing", con)
  close(con)
}

# ---- Nested failure to actually crash ----------------------------------------
third <- function(msg, code) stop(msg, call. = TRUE)
second <- function(msg, code) third(msg, code)
first <- function(msg, code) second(msg, code)

first(msg, as.integer(ecode))
