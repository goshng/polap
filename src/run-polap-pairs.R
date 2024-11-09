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

# name: select long reads
#
# synopsis:
#	run-polap-pairs.R mt.contig.name-1 contig.tab ${MTSEEDSDIR} $SINGLE_MIN $SINGLE_MIN
#
# requirement: executes Flye 
# flye --nano-raw "$LR3K" --out-dir "$ODIR" \
# 	--threads "${_arg_threads}" \
# 	--stop-after contigger \
# 	--asm-coverage 30 \
# 	--genome-size "$EXPECTED_GENOME_SIZE"
#
# input: 
#   1. mt.contig.name-1 - contig or edge numbers
#   2. contig.tab - minimap2 output modified
#   3. seeds directory for output
#   4. pair minimum length V11: 3000 for MT, 1000 for PT
#   5. bridge minimum length V7: depends: 3000 or 5000
#   6. single minimum length V11: 3000 for MT, 0 or 1000 for PT
#
# output: 
#   single.names and <contig1-contig2.name> files in the output directory
#
# MTSEEDSDIR="$ODIR"/60-mt-${STEP4}/o${MR}/seeds
#	run-polap-pairs.R mt.contig.name-1 contig.tab ${MTSEEDSDIR} $PAIR_MIN $BRIDGE_MIN $SINGLE_MIN
#	run-polap-pairs.R mt.contig.name-1 contig.tab ${MTSEEDSDIR} $PAIR_MIN $BRIDGE_MIN

suppressPackageStartupMessages(library("dplyr"))

# see https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/
# see https://blog.sellorm.com/2017/12/18/learn-to-write-command-line-utilities-in-r/
args = commandArgs(trailingOnly=TRUE)

x = as_tibble(read.table(args[1]))
if (nrow(x) > 1) {

y = t(combn(x$V1,2))

x = as_tibble(read.table(args[2]))
mtdir = args[3]
pair_min = as.numeric(args[4])
brigde_min = as.numeric(args[5])
#single_min = as.numeric(args[5])
if (length(args) == 6) {
  single_min = as.numeric(args[6])
} else {
  single_min = 0
}
stopifnot(pair_min >= 0, brigde_min >= 0, single_min >= 0)

# https://lh3.github.io/minimap2/minimap2.html
# Li 2008 - "OUTPUT FORMAT: Minimap2 outputs mapping positions in the Pairwise mApping Format (PAF) by default. PAF is a TAB-delimited text format with each line consisting of at least 12 fields as are described in the following table:"
orders <- function(y){
  stopifnot(length(y) == 2)
  z = as.numeric(gsub("\\D", "", y))
  bname = paste0(z[1],"_",z[2])
  oname = paste0(mtdir,"/",bname,".name")

  y1 = x %>% filter(V6==y[1], V10/V11 > 0.7, V11 > pair_min, V7 > brigde_min) %>% select(V1) %>% distinct(V1)
  y2 = x %>% filter(V6==y[2], V10/V11 > 0.7, V11 > pair_min, V7 > brigde_min) %>% select(V1) %>% distinct(V1)

  intersect(y1,y2) %>% write.table(oname,row.names=F,col.names=F,quote=F)

}

sfile = paste0(mtdir,"/single.names")
# x %>% filter((V4-V3)/V2 > 0.7, V10/V11 > 0.7) %>% select(V1) %>% distinct(V1) %>% write.table(sfile,row.names=F,col.names=F,quote=F)
x %>% filter((V4-V3)/V2 > 0.7, V10/V11 > 0.7, V11 > single_min) %>% select(V1) %>% distinct(V1) %>% write.table(sfile,row.names=F,col.names=F,quote=F)

apply(y, 1, orders)

} else {

x = as_tibble(read.table(args[2]))
mtdir = args[3]
single_min = as.numeric(args[4])

sfile = paste0(mtdir,"/single.names")
# x %>% filter((V4-V3)/V2 > 0.7, V10/V11 > 0.7) %>% select(V1) %>% distinct(V1) %>% write.table(sfile,row.names=F,col.names=F,quote=F)
x %>% filter(V10/V11 > 0.7, V11 > single_min) %>% select(V1) %>% distinct(V1) %>% write.table(sfile,row.names=F,col.names=F,quote=F)

}
