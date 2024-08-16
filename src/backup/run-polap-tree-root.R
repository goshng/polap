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

library(ggtree)
library(treeio)
library(readr)

tree = read.newick("phylogenies/s.treefile", node.label = "support")
xm = read_tsv("species.map", col_names = F)

tree4 = root(tree,2)

tree2 = tree_subset(tree4,48,levels_back = 0)
root <- rootnode(tree2)
tree3 = rename_taxa(tree2, xm, X1, X2)

ggtree(tree, branch.length="none") + 
  geom_text(aes(label=node), hjust=-.3) +
  geom_tiplab(offset=0.5)

#  xlim(-0.01, 0.35)    

ggtree(tree3, color="black", size=0.5, linetype=1,  right=T, layout="roundrect") +
    geom_text2(aes(subset=!isTip, label=node), size=5.1) +
  geom_rootedge(rootedge = 0.05) +
  geom_tiplab(size=3, hjust = -0.060, fontface="italic", align = T) +
  hexpand(.35) +
  geom_treescale(x=0.01, y=0.5, fontsize=3, linesize=0.5, offset=0.5) +
  geom_point2(aes(subset=!isTip & node != root, 
                  fill=cut(support, c(0, 90, 95, 100))), 
              shape=21, size=3) + 
  scale_fill_manual(values=c("white", "grey", "black"), guide='legend', 
                    name='Bootstrap value', 
                    breaks=c('(95,100]', '(90,95]', '(0,90]'), 
                    labels=expression(BP>95, 90 < BP * " <= 95", BP <= 90))



  #xlim(-0.01, 0.35) +  

# backup

g = ggtree(tree3, color="black", size=0.5, linetype=1,  right=T, layout="roundrect") +
  geom_text2(aes(subset=!isTip, label=node), size=0.1) +
  geom_rootedge(rootedge = 0.01) +
  geom_tiplab(size=3, hjust = -0.060, fontface="italic", align = T) +
  hexpand(.35) +
  geom_treescale(x=0.01, y=0.5, fontsize=3, linesize=0.5, offset=0.5) +
  geom_point2(aes(subset=!isTip & node != root, 
                  fill=cut(support, c(0, 90, 95, 100))), 
              shape=21, size=3) + 
  scale_fill_manual(values=c("white", "grey", "black"), guide='legend', 
                    name='Bootstrap value', 
                    breaks=c('(95,100]', '(90,95]', '(0,90]'), 
                    labels=expression(BP>95, 90 < BP * " <= 95", BP <= 90))



pdf("t2.pdf", width = 11, height = 7)
g
dev.off()

print("Check: open t2.pdf")