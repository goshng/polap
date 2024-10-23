#!/usr/bin/env Rscript
library("optparse")
parser <- OptionParser()
parser <- add_option(parser, c("-v", "--verbose"), action="store_true",
                default=TRUE, help="Print extra output [default]")
parser <- add_option(parser, c("-q", "--quietly"), action="store_false",
                    dest="verbose", help="Print little output")
parser <- add_option(parser, c("-c", "--count"), type="integer", default=5,
                help="Number of random normals to generate [default %default]",
                metavar="number")
parse_args(parser) #, args = c("--quietly", "--count=15"))



