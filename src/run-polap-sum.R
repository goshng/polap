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

flye_dir = "o"
if (exists("out")) {
    flye_dir = out
}


assembly_info = paste0(flye_dir,"/assembly_info.txt")
mt.gene.count = paste0(flye_dir,"/50-annotation/mt.gene.count")
pt.gene.count = paste0(flye_dir,"/50-annotation/pt.gene.count")

if (!file.exists(assembly_info) || 
    !file.exists(mt.gene.count) ||
    !file.exists(pt.gene.count)) {
    print("ERROR: some files:", assembly_info, mt.gene.count, pt.gene.count, "does not exist!")
}

x = as_tibble(read.table(assembly_info))
y.mt = as_tibble(read.table(mt.gene.count)) %>% rename(mt=V2)
y.pt = as_tibble(read.table(pt.gene.count)) %>% rename(pt=V2)
z = x %>% left_join(y.mt) %>% left_join(y.pt) %>% arrange(desc(mt)) %>% filter(!is.na(mt)) %>% filter(mt > 0 | pt > 0)
z.1 = x %>% left_join(y.mt) %>% left_join(y.pt) 
z.v = x %>% left_join(y.mt) %>% left_join(y.pt) %>% arrange(mt) %>% filter(is.na(mt))
print("print variable: z")
