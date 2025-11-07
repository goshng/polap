# Version: v0.1.0
options(error = function(e) {
  n <- sys.nframe()
  call <- if (n > 0) deparse(sys.call(n)) else "?"
  msg <- conditionMessage(e)
  cat(sprintf("R-FAILED ?:?: %s\n", msg), file = stderr())
  utils::traceback() # full stack
  q(status = 10, save = "no")
})
