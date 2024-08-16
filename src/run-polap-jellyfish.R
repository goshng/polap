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

args = commandArgs(trailingOnly=TRUE)

x <- read.table(args[1])
total_read_length = as.numeric(args[2])

n <- nrow(x)
a <- which.min(x$V2[1:10])
b <- n - 10 + which.min(x$V2[(n-9):n])
y <- which.max(x$V2[a:b]) + (a - 1)
z <- sum(as.numeric(x[a:b,1] * x[a:b,2]))/y

write.table(total_read_length/z, args[3],row.names=F,col.names=F,quote=F)
write.table(z, args[4],row.names=F,col.names=F,quote=F)

