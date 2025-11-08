# polaplib/polap_r_failfast.R
# Make R failures print a readable stack then exit with a code (default 7).

options(keep.source = TRUE, warn = 1, show.error.messages = FALSE)

polap_r_print_stack <- function() {
  calls <- sys.calls()
  if (!length(calls)) { cat("Stack (empty)\n"); return(invisible()) }
  cat("Stack (oldest → newest):\n")
  for (i in seq_along(calls)) {
    cl <- calls[[i]]
    sr <- attr(cl, "srcref")
    txt <- paste(deparse(cl), collapse = " ")
    if (!is.null(sr)) {
      sf <- attr(sr, "srcfile")
      fn <- if (!is.null(sf) && !is.null(sf$filename)) basename(sf$filename) else "<src>"
      cat(sprintf("#%d %s:%d:%d  %s\n", i, fn, as.integer(sr[[1]]), as.integer(sr[[2]]), txt))
    } else {
      cat(sprintf("#%d %s\n", i, txt))
    }
  }
  invisible()
}

polap_r_enable_failsafe <- function(exit_code = 7L, with_stack = TRUE) {
  if (isTRUE(with_stack)) {
    options(error = function() {
      try(polap_r_print_stack(), silent = TRUE)
      quit(save = "no", status = as.integer(exit_code))
    })
  } else {
    options(error = function() quit(save = "no", status = as.integer(exit_code)))
  }
}

