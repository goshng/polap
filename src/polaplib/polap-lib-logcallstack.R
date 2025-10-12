#!/usr/bin/env Rscript
# Version: v0.3.5 (caller is two frames up; robust file/line)

if (!isTRUE(getOption("keep.source"))) options(keep.source = TRUE)

.polap_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
.polap_script <- (function() {
  ca <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", ca[grep("^--file=", ca)])
  if (length(f) > 0) basename(f) else "Rscript"
})()

# Safe helpers for filename/line
.polap_fn_file <- function(fn) {
  val <- tryCatch(utils::getSrcFilename(fn, full.names = FALSE),
    error = function(e) character(0)
  )
  if (length(val) == 0L || isTRUE(is.na(val))) .polap_script else basename(val)
}
.polap_fn_line <- function(fn) {
  ln <- tryCatch(utils::getSrcLocation(fn, first = TRUE),
    error = function(e) NA_integer_
  )
  if (length(ln) == 0L) NA_integer_ else ln
}
.polap_call_line <- function(i) {
  sc <- sys.calls()
  if (i >= 1L && i <= length(sc)) {
    sr <- attr(sc[[i]], "srcref")
    if (!is.null(sr)) {
      ln <- sr[[1]]
      if (length(ln) > 0L) {
        return(ln)
      }
    }
  }
  fn <- tryCatch(sys.function(i), error = function(e) NULL)
  if (is.null(fn)) NA_integer_ else .polap_fn_line(fn)
}

# ---------- FIX: header uses the *caller of the caller* (sys.parent(2)) ----------
.polap_header <- function() {
  i <- sys.parent(2L) # caller of polap_log_* (skip logger frame)
  if (i == 0L) i <- sys.parent(1L) # fallback: immediate parent
  if (i == 0L) i <- sys.nframe() # fallback: top
  fn_name <- tryCatch(as.character(sys.call(i)[[1]]), error = function(e) "main")
  fn_obj <- tryCatch(sys.function(i), error = function(e) NULL)
  file <- if (!is.null(fn_obj)) .polap_fn_file(fn_obj) else .polap_script
  line <- .polap_call_line(i)
  sprintf("[%s %s@%s:%s]", .polap_now(), fn_name, file, ifelse(is.na(line), "NA", line))
}

# Stack printer (skip logger internals)
polap_stack_dump <- function() {
  n <- sys.nframe()
  if (n <= 0L) {
    return(invisible(NULL))
  }
  skip <- c("polap_stack_dump", ".polap_header", "polap_log_err", "polap_log_warn", "polap_log_info")
  k <- 0L
  for (i in rev(seq_len(n))) {
    fn_name <- tryCatch(as.character(sys.call(i)[[1]]), error = function(e) "main")
    if (fn_name %in% skip) next
    fn_obj <- tryCatch(sys.function(i), error = function(e) NULL)
    file <- if (!is.null(fn_obj)) .polap_fn_file(fn_obj) else .polap_script
    line <- .polap_call_line(i)
    cat(
      sprintf(
        "[%s %s@%s:%s] #%d\n",
        .polap_now(), fn_name, file, ifelse(is.na(line), "NA", line), k
      ),
      file = stderr()
    )
    k <- k + 1L
  }
}

# Logging API
polap_log_info <- function(...) {
  hdr <- .polap_header()
  cat(hdr, paste(...), "\n", sep = " ")
}
polap_log_warn <- function(...) {
  hdr <- .polap_header()
  cat(hdr, "WARNING:", paste(...), "\n", sep = " ", file = stderr())
}
polap_log_err <- function(...) {
  hdr <- .polap_header()
  cat(hdr, "ERROR:", paste(...), "\n", sep = " ", file = stderr())
  polap_stack_dump()
}

# Optional: single top-level trap (don’t mix with tryCatch around the same call)
polap_enable_top_error_trap <- function(status = 1L) {
  options(error = function(e) {
    polap_log_err(conditionMessage(e))
    quit(save = "no", status = status)
  })
}

# Optional: source with keep.source=TRUE
polap_source <- function(file, local = TRUE) {
  base::source(file, local = local, keep.source = TRUE, echo = FALSE, chdir = FALSE)
}
