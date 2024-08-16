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

tree = read.newick("phylogenies/s2.treefile", node.label = "support")

pdf("t1.pdf", width = 11, height = 7)
ggtree(tree, branch.length="none") + 
  geom_text(aes(label=node), hjust=-.3) +
  geom_tiplab(offset=0.5)
dev.off()

print("Check: open t1.pdf")