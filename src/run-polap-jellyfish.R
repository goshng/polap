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

# name: estimates genome size using Jellyfish output
#
# synopsis:
# run-polap-jellyfish.R <jellyfish .histo file> <long read total length> <out:long_coverage> <out:short_expected_genome_size>
# run-polap-jellyfish.R 19mer_out.histo 2092086846 long_coverage.txt short_expected_genome_size.txt
#
# requirement: executes Jellyfish 
# jellyfish count -t 4 -C -m 19 -s 5G -o 19mer_out --min-qual-char=? s1.fq s2.fq
# jellyfish histo -o 19mer_out.histo 19mer_out
#
# input: Jellyfish .histo output file
# output: coverage and genome size
#   "$WDIR"/run-polap-jellyfish.R 19mer_out.histo "$LONG_TOTAL_LENGTH" coverage.txt expected_genome_size.txt
#	EXPECTED_GENOME_SIZE=$(cat expected_genome_size.txt)
#	EXPECTED_GENOME_SIZE=${EXPECTED_GENOME_SIZE%.*}
#	EXPECTED_COVERAGE=$(cat coverage.txt)
#	EXPECTED_COVERAGE=${EXPECTED_COVERAGE%.*}
#	echo "DATA: short reads expected genome size (bases): $EXPECTED_GENOME_SIZE"
#	echo "DATA: long reads expected coverage: ${EXPECTED_COVERAGE}x"

args <- commandArgs(trailingOnly = TRUE)

# two input arguments
x <- read.table(args[1])
total_read_length <- as.numeric(args[2])

n <- nrow(x)
a <- which.min(x$V2[1:10])
b <- n - 10 + which.min(x$V2[(n - 9):n])
y <- which.max(x$V2[a:b]) + (a - 1)
z <- sum(as.numeric(x[a:b, 1] * x[a:b, 2])) / y

# two output argments
write.table(total_read_length / z, args[3],
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(z, args[4], row.names = FALSE, col.names = FALSE, quote = FALSE)