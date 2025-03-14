#!/usr/bin/env Rscript
################################################################################
# This file is part of polap.
#
# polap is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# polap is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with
# polap. If not, see <https://www.gnu.org/licenses/>.
################################################################################

args <- commandArgs(trailingOnly = TRUE)

# Test if the number of arguments is exactly 2
if (length(args) != 2) {
  # stop("Error: Exactly two command line arguments are required.")
  # args <- c("Pisum_sativum/jellyfish_out.histo", "Pisum_sativum/short_expected_genome_size.txt")
}

tryCatch(
  {
    # input arguments
    x <- read.table(args[1])
    n <- nrow(x)
    a <- which.min(x$V2[1:10])
    b <- n - 10 + which.min(x$V2[(n - 9):n])
    y <- which.max(x$V2[a:b]) + (a - 1)
    z <- sum(as.numeric(x[a:b, 1] * x[a:b, 2])) / y

    # output argments
    # The maximum integer (.Machine$integer.max) is 2,147,483,647.
    # more than 2Gb genome size would print an error.
    # write.table(as.integer(z), args[2], row.names = FALSE, col.names = FALSE, quote = FALSE)
    # So, we use the following:
    write(format(trunc(z), scientific = FALSE, trim = TRUE), file = args[2])
  },
  error = function(e) {
    # cat("Error: ", args[0], ": ", e$message, "\n")
    quit(status = 1) # Exit with status 1 (indicating an error)
  }
)
