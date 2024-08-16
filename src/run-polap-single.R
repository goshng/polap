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

mtdir = args[2]
single_min = as.numeric(args[3])

sfile = paste0(mtdir,"/single.names")
x %>% filter(V10/V11 > 0.7,V11 > single_min) %>% select(V1) %>% distinct(V1) %>% write.table(sfile,row.names=F,col.names=F,quote=F)

