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

suppressPackageStartupMessages(library("dplyr"))
args = commandArgs(trailingOnly=TRUE)


x = as_tibble(read.table(args[1]))
copy_number_min = median(x$V6)
# x %>% filter(V2<3000000,V2>1000) %>% filter(V6>copy_number_min) %>% select(V1) %>% write.table(args[2],row.names=F,col.names=F,quote=F)
x %>% filter(V2>1000,V6>copy_number_min) %>% select(V1) %>% write.table(args[2],row.names=F,col.names=F,quote=F)
