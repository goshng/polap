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

# name: phylogeny and gene presence/absence table of Carex pseudochinensis mtDNA
#
# synopsis:
#
# requirement:
#   amino acid sequence files
#   orthofinder
#   iqtree
#
# input:
#   s.treefile
#   species.map
#   genes.tsv
#   outgroup
#   outgroup node
#   root node
#
# output:
#   figure

suppressPackageStartupMessages(library("ggtree"))
suppressPackageStartupMessages(library("treeio"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("aplot"))
suppressPackageStartupMessages(library("stringr"))

#opt <- list("tree" = "s1.treefile", "node" = 2, "root" = 10, "map" = "species.map", "genes" = "genes.tsv", "outgroup" = "NC_037304.1") 
opt <- list("tree" = "s.treefile", "map" = "species.map", "genes" = "genes.tsv", "outgroup" = "NC_037304") 
opt <- list("tree" = "phylogenies/s1.treefile", "map" = "species.map", "genes" = "phylogenies/genes.tsv", "outgroup" = "NC_037304") 

# --tree phylogenies/s2.treefile
if (!is.null(opt$tree)) {
  tree = read.newick(opt$tree, node.label = "support")
}

# --map species.map --genes phylogenies/genes.csv --outgroup NC_036467,NC_044153,NC_073519
if (!is.null(opt$map)) {
  xm = read_tsv(opt$map, col_names = F, show_col_types = FALSE)
}

# --map species.map --genes phylogenies/genes.csv --outgroup NC_036467,NC_044153,NC_073519
if (!is.null(opt$genes) && !is.null(opt$outgroup)) {
  x = strsplit(opt$outgroup,split=",",fixed=TRUE)[[1]]
  xl = read_tsv(opt$genes, show_col_types = FALSE) %>% select(!all_of(x))
}

tree %>% print(n=50)

if (is.null(opt$node)) {
  tree2 = tree        
} else {
  tree2 = root(tree,opt$node)
}

tree2 %>% print(n=50)

ggtree(tree2) +
  geom_text(aes(label=node), hjust=-.3) +
  geom_tiplab(offset = 0.005)

ggtree(tree2, branch.length="none") + 
  geom_text(aes(label=node), hjust=-.3) +
  geom_tiplab(offset=0.5)


if (is.null(opt$node)) {
  tree4 = tree        
} else {
  tree4 = root(tree,outgroup=opt$node)
}

tree4 %>% print(n=50)


#tree2 = tree_subset(tree4, opt$node, levels_back = 5)
#tree2 %>% print(n=20)
#tree2 = tree_subset(tree4, opt$root, levels_back = 0)
root <- rootnode(tree2)

tree3 = rename_taxa(tree2, xm, X1, X2)
rootnode(tree3)
tree3 %>% print(n=20)

if (is.null(opt$nodesize)) {
  nodesize = 0.1
} else {
  nodesize = 5.1
}

# S. cypria
# xlim(-0.01, 0.35) +  

# hexpand(.35) +

# Carex pseudochinensis
# xlim(-0.01, 1.3) +  # CHANGE
ggtree(tree3, color="black", size=0.5, linetype=1,  right=T, layout="roundrect") +
  geom_text2(aes(subset=!isTip, label=node), size=nodesize) +
  geom_rootedge(rootedge = 0.01)
  
g = ggtree(tree3, color="black", size=0.5, linetype=1,  right=T, layout="roundrect") +
  geom_text2(aes(subset=!isTip, label=node), size=nodesize) +
  geom_rootedge(rootedge = 0.01) +
  geom_tiplab(size=3, hjust = -0.060, fontface="italic", align = T) +
  hexpand(.35) +
  xlim(-0.01, 0.7) +  # CHANGE
  geom_treescale(x=0.01, y=0.5, fontsize=3, linesize=0.5, offset=0.5) +
  geom_point2(aes(subset=!isTip & node != root, 
                  fill=cut(support, c(0, 90, 95, 100))), 
              shape=21, size=3) + 
  scale_fill_manual(values=c("white", "grey", "black"), guide='legend', 
                    name='Bootstrap value', 
                    breaks=c('(95,100]', '(90,95]', '(0,90]'), 
                    labels=expression(BP>95, 90 < BP * " <= 95", BP <= 90))
g1 = g

g1

# step3

g2 = g1 + 
  geom_cladelab(node=12, label="Cyperaceae", align=TRUE, offset = 0.35, offset.text = 0.005) +
  geom_cladelab(node=17, label="Juncaceae", align=TRUE, offset = 0.35, offset.text=0.005) 
g2

#g2
# step4
# xm = read_tsv("species.map", col_names = F)
# xl = read_csv("phylogenies/genes.csv") %>% select(!c("NC_036467","NC_044153","NC_073519"))

# "chr1_cds-blatx_"
#target_gene = "_cds-blatx_"
target_gene = "cpseudochinensis_"

# tl = gsub(target_gene,"", xl$target)
# tl = gsub("_1","", tl)

last.target = xl[,ncol(xl)]
tl = str_extract(last.target[[1]], paste0(target_gene, "([a-zA-Z0-9]+)"), group=1)
tl[is.na(tl)] = xl[is.na(tl),ncol(xl)][[1]]

last.target[[1]]

# single target
#target_gene = paste0(opt$locus, "_cds-blatx_")
#tl = gsub(target_gene,"", xl$target)
#tl = gsub("_1","", tl)

ngene = nrow(xl)
nspecies = ncol(xl)
xl2 = xl
xl[!is.na(xl)] <- "O"
xl[is.na(xl)] <- "X"

dl = data.frame(label=rep(colnames(xl),each=ngene), 
                category=rep(tl,nspecies), 
                value=as.numeric(xl=="O"))

#grepl(target_gene, last.target[[1]])

pos.lack.target = c(rep(F,(nspecies-1)*ngene), !grepl(target_gene, last.target[[1]]))
dl$value[pos.lack.target] = 0
dl$value[grepl(",",as.matrix(xl2))] = 2 # no treat for three or more commas 

dl2 = dl %>% mutate(label2=xm$X2[match(unlist(dl$label), xm$X1)]) %>% select(!label) %>% rename(label=label2)

p2 <- ggplot(dl2, aes(x=category, y=label, fill=value)) + 
  geom_tile(color="black", show.legend = F) + 
  scale_fill_gradient(low="white", high="grey") + 
  theme(axis.text.x = element_text(angle=90,size=6,hjust = 0.1,vjust=0.2)) + 
  scale_y_discrete(labels = NULL, breaks = NULL) + labs(y = "") +
  xlab(NULL) + ylab(NULL)

p2


# step5

toffset <- 0.35
g2 = g1 + 
  geom_cladelab(node=12, label="Cyperaceae", align=TRUE, offset = toffset, offset.text = 0.005) +
  geom_cladelab(node=17, label="Juncaceae", align=TRUE, offset = toffset, offset.text=0.005) 



pdf("figure3.pdf", width = 11, height = 7)
g <- p2 %>% insert_left(g2, width=1.2)
print(g)
dev.off()

cat("Check: open", "figure3.pdf", "\n")

#g2

