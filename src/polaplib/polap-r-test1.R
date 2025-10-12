#!/usr/bin/env Rscript
# polap-r-test1.R
base <- Sys.getenv("POLAPLIB_DIR", Sys.getenv("_POLAPLIB_DIR", "."))
source(file.path(base, "polap-lib-logcallstack.R"), local = TRUE)
# EITHER enable a global trap...
polap_enable_top_error_trap(status = 2)

inner_r <- function() {
  polap_log_err("R: simulated failure inside inner_r()")
  stop("boom from R")
}
outer_r <- function() {
  polap_log_info("R: entering outer_r()")
  inner_r()
}

# ...and just call outer_r()
outer_r()

# If you prefer tryCatch instead of a global trap, then:
# tryCatch(outer_r(),
#   error = function(e) {
#     polap_log_err(paste("R: unhandled error at top:", conditionMessage(e)))
#     quit(save="no", status=2)
#   }
# )
