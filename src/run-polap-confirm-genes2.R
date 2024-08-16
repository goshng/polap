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
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("stringr"))

# args[1] = "genes.csv"
# args[2] = "genes.target.txt"

x = read_tsv(args[1], show_col_types = FALSE)

# y = unlist(unique(as.vector(x[35,])))
# z = y %>% str_split(",") %>% unlist() %>% unique() %>% na.omit()
# z = gsub("-.+","",z) 
# z = gsub("\\[gene=","", z)
# z = gsub("\\]","", z) %>% unique()

# paste(z, collapse = ",") 

# gsub("-.+","",y) %>% gsub("gene=","") %>% str_split(",") %>% unlist() %>% unique() %>% na.omit()

u1 <- function(x) {
  y = unlist(unique(as.vector(x)))
  z = y %>% str_split(",") %>% unlist() %>% unique() %>% na.omit()
  z = gsub("-.+","",z) 
  z = gsub("\\[gene=","", z)
  z = gsub("\\]","", z) %>% unique()
  z[1]
#  paste(z, collapse = ",") 
}

x1 = apply(x,1,u1)

g = read_tsv(args[2],skip_empty_rows = F,show_col_types = FALSE)
g[is.na(g),] = x1[is.na(g)]

cbind(x,g) %>% write.table(args[3],row.names=F,col.names=T,quote=F,sep="\t")

