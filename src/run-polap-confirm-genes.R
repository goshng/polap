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
suppressPackageStartupMessages(library("tidyr"))

# id = gsub(".txt","",basename(args[3]))
id = args[3]

x = as_tibble(read.table(args[1],header=F,sep="\t",skip = 1)) %>% rename(id=V2) %>% separate_longer_delim(c(id), delim = ", ")
y = as_tibble(read.table(args[2])) %>% rename(id=V1)
# z = left_join(x,y,by='id') %>% group_by(V1) %>% reframe(V2 = paste(V2, collapse = ",")) %>% select(V2) %>% rename(!!id := V2)
z = left_join(x,y,by='id') %>% group_by(V1) %>% reframe(V2 = paste(V2, collapse = ",")) %>% rename(!!id := V2)


z %>% write.table(args[4],row.names=F,col.names=T,quote=F,sep="\t")
