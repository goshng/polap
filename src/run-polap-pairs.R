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

suppressPackageStartupMessages(library("dplyr"))

x = as_tibble(read.table(args[1]))
y = t(combn(x$V1,2))

x = as_tibble(read.table(args[2]))
mtdir = args[3]
pair_min = as.numeric(args[4])
single_min = as.numeric(args[5])

orders <- function(y){
  stopifnot(length(y) == 2)
  z = as.numeric(gsub("\\D", "", y))
  bname = paste0(z[1],"_",z[2])
  oname = paste0(mtdir,"/",bname,".name")

# v1
  # y1 = x %>% filter(V6==y[1]) %>% select(V1) %>% distinct(V1)
  # y2 = x %>% filter(V6==y[2]) %>% select(V1) %>% distinct(V1)
# v3
  y1 = x %>% filter(V6==y[1], V10/V11 > 0.7, V11 > pair_min) %>% select(V1) %>% distinct(V1)
  y2 = x %>% filter(V6==y[2], V10/V11 > 0.7, V11 > pair_min) %>% select(V1) %>% distinct(V1)
# v2
  # y1 = x %>% filter(V10/V11 > 0.7) %>% filter(V6==y[1]) %>% select(V1) %>% distinct(V1)
  # y2 = x %>% filter(V10/V11 > 0.7) %>% filter(V6==y[2]) %>% select(V1) %>% distinct(V1)

  intersect(y1,y2) %>% write.table(oname,row.names=F,col.names=F,quote=F)

}

# version 0.1.1: 0.9 not 0.7
sfile = paste0(mtdir,"/single.names")
# v1
# x %>% filter((V4-V3)/V2 > 0.9) %>% select(V1) %>% distinct(V1) %>% write.table(sfile,row.names=F,col.names=F,quote=F)
# v2
# x %>% filter((V4-V3)/V2 > 0.7) %>% select(V1) %>% distinct(V1) %>% write.table(sfile,row.names=F,col.names=F,quote=F)
# x %>% filter((V4-V3)/V2 > 0.7,V11 >= 1000) %>% select(V1) %>% distinct(V1) %>% write.table(sfile,row.names=F,col.names=F,quote=F)
x %>% filter((V4-V3)/V2 > 0.7, V10/V11 > 0.7, V11 > single_min) %>% select(V1) %>% distinct(V1) %>% write.table(sfile,row.names=F,col.names=F,quote=F)

apply(y, 1, orders)

# Carex siderosticta mtDNA case:
# (m1) [goshng@Gingko cpdna]$ PAIR_MIN=10000 SINGLE_MIN=10000 POLAP_SKIP=1 polap -4 2 -k n3k.fq.gz
# 5000 trying
# this value could be or should be greater than 3 to 4 times of gene length (1kb)

#  1075  PAIR_MIN=9000 SINGLE_MIN=7000 POLAP_SKIP=1 polap -4 4 -k n3k.fq.gz &
#  1076  PAIR_MIN=7000 SINGLE_MIN=9000 POLAP_SKIP=1 polap -4 5 -k n3k.fq.gz &


#  1073  PAIR_MIN=7000 SINGLE_MIN=7000 POLAP_SKIP=1 polap -4 3 -k n3k.fq.gz &
#  1052  PAIR_MIN=5000 SINGLE_MIN=5000 POLAP_SKIP=1 polap -4 2 -k n3k.fq.gz
#  1011  PAIR_MIN=4000 SINGLE_MIN=4000 POLAP_SKIP=1 polap -4 8 -k n3k.fq.gz &
#  1013  PAIR_MIN=3000 SINGLE_MIN=6000 POLAP_SKIP=1 polap -4 9 -k n3k.fq.gz &
#  1010  PAIR_MIN=3000 SINGLE_MIN=3000 POLAP_SKIP=1 polap -4 7 -k n3k.fq.gz &
#  1009  PAIR_MIN=2000 SINGLE_MIN=2000 POLAP_SKIP=1 polap -4 6 -k n3k.fq.gz &

# c. sider... - 5000 is okay not 4000 or below
# 3000/6000?