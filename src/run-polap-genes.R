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

#!/var2/home/alice/miniconda3/envs/bio38/bin/Rscript --vanilla

suppressPackageStartupMessages(library("dplyr"))


x = as_tibble(read.table(args[1]))
y1 = x %>% select(V2,V9,V10) %>% filter(V9<V10) %>% rename(V1=V2,V2=V9,V3=V10)
y2 = x %>% select(V2,V9,V10) %>% filter(V9>V10) %>% rename(V1=V2,V2=V10,V3=V9)

#y1 = x %>% filter(V2=="contig_19") %>% select(V2,V9,V10) %>% filter(V9<V10) %>% rename(V1=V2,V2=V9,V3=V10)
#y2 = x %>% filter(V2=="contig_19") %>% select(V2,V9,V10) %>% filter(V9>V10) %>% rename(V1=V2,V2=V10,V3=V9)

z = bind_rows(y1,y2)

write.table(z, args[2],row.names=F,col.names=F,quote=F,sep='\t')

