#!/usr/bin/env Rscript
# Version: v0.1.0
f <- function() g()
g <- function() stop("demo crash: R") # prints via our R profile handler
f()
