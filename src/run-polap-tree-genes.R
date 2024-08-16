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

library(ggplot2)
library(readr)
library(dplyr)

s = "AJ506156,AP000423,GU195652"
x = strsplit(s,split=",",fixed=TRUE)[[1]]
xl = read_tsv("phylogenies/genes.tsv", show_col_types = FALSE) %>% select(!all_of(x))


xm = read_tsv("species.map", col_names = F)
xl = read_csv("phylogenies/genes.csv") %>% select(!c("NC_036467","NC_044153","NC_073519"))

tl = gsub("chr1_cds-blatx_","", xl[,ncol(xl)])

tl = gsub("_1","", tl)
ngene = nrow(xl)
nspecies = ncol(xl)
xl2 = xl
xl[!is.na(xl)] <- "O"
xl[is.na(xl)] <- "X"

dl = data.frame(label=rep(colnames(xl),each=ngene), 
                category=rep(tl,nspecies), 
                value=as.numeric(xl=="O"))

pos.lack.target = c(rep(F,(nspecies-1)*ngene), !grepl("chr1_cds-blatx_", xl2$target))
dl$value[pos.lack.target] = 0
dl$value[grepl(",",as.matrix(xl2))] = 2 # no treat for three or more commas 


dl2 = dl %>% mutate(label2=xm$X2[match(unlist(dl$label), xm$X1)]) %>% select(!label) %>% rename(label=label2)

p2 <- ggplot(dl2, aes(x=category, y=label, fill=value)) + 
  geom_tile(color="black", show.legend = F) + 
  scale_fill_gradient(low="white", high="grey") + 
  theme(axis.text.x = element_text(angle=90,size=6,hjust = 0.1,vjust=0.2)) + 
  scale_y_discrete(labels = NULL, breaks = NULL) + labs(y = "") +
  xlab(NULL) + ylab(NULL)



pdf("p2.pdf", width = 11, height = 7)
p2
dev.off()

print("Check: open p2.pdf")

