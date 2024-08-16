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

# see: CHANGE
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stats"))

option_list <- list( 
    make_option("--tree", help = "Newick tree"),
    make_option("--pdf", default="o.pdf", help = "PDF output file [default \"%default\"]"),
    make_option("--map", help = "species.map of taxon to species"),
    make_option("--genes", help = "genes.csv"),
    make_option("--locus", help = "locus tag"),
    make_option("--outgroup", help = "outgroups: NC_036467,NC_044153,NC_073519"),
    make_option("--node", type = "integer", help = "node number for outgroup"),
    make_option("--nodesize", action="store_false", help = "label size"),
    make_option("--root", type = "integer", help = "root node number"),
    make_option("--step", default=0, help = "Step [default \"%default\"]")
    )
                                        
opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library("ggtree"))
suppressPackageStartupMessages(library("treeio"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("aplot"))
suppressPackageStartupMessages(library("stringr"))


step1 <- function() {

    tree = read.newick("phylogenies/s2.treefile", node.label = "support")
    # print(tree)

    if (is.null(opt$node)) {
        tree2 = tree        
    } else {
        tree2 = root(tree,opt$node)
    }


# Note: http://www.sthda.com/english/wiki/ggsave-save-a-ggplot-r-software-and-data-visualization
    pdf(opt$pdf, width = 11, height = 7)
    # pdf("t1.pdf", width = 11, height = 7)
    g = ggtree(tree2, branch.length="none") + 
        geom_text(aes(label=node), hjust=-.3) +
        geom_tiplab(offset=0.5)
    print(g)
    dev.off()

    cat("Check: open", opt$pdf, "\n")
}


step2 <- function(draw=TRUE) {

    # tree = read.newick("phylogenies/s2.treefile", node.label = "support")
    # xm = read_tsv("species.map", col_names = F)
    # rootnodenumber = 26

    if (is.null(opt$node)) {
        tree4 = tree        
    } else {
        tree4 = root(tree,opt$node)
    }

    tree2 = tree_subset(tree4, opt$root, levels_back = 0)
    root <- rootnode(tree2)
    tree3 = rename_taxa(tree2, xm, X1, X2)

    

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

    if (draw) {
        pdf(opt$pdf, width = 11, height = 7)
        print(g)
        dev.off()
        cat("Check: open", opt$pdf, "\n")
    }
    return(g)
}

step5_scutellaria_cypria <- function(g,p2,draw=TRUE) {

    # tree = read.newick("phylogenies/s2.treefile", node.label = "support")
    # xm = read_tsv("species.map", col_names = F)
    # rootnodenumber = 26


    g2 = g %>% flip(27,31) %>% rotate(29) +
        geom_cladelab(node=36, size=1, label="Oleaceae", align=TRUE, offset = 0.15, offset.text=0.005) +
        geom_cladelab(node=14, label="Gesneriaceae", align=TRUE, offset = 0.15, offset.text=0.005) +
        geom_cladelab(node=25, label="Plantaginaceae", align=TRUE, offset = 0.15, offset.text=0.005) +
        geom_cladelab(node=13, label="Phrymaceae", align=TRUE, offset = 0.15, offset.text=0.005) +
        geom_cladelab(node=32, label="Orobanchaceae", align=TRUE, offset = 0.15, offset.text=0.005) +
        geom_cladelab(node=27, label="Lamiaceae", align=TRUE, offset = 0.15, offset.text=0.005)


    pdf(opt$pdf, width = 11, height = 7)
    p2 %>% insert_left(g2, width=1.8) %>% print()
    # print(g)
    dev.off()
    cat("Check: open", opt$pdf, "\n")
}

# C. pseudochinensis
step3_carex_pseudochinensis <- function(g1,draw=TRUE) {

    # tree = read.newick("phylogenies/s2.treefile", node.label = "support")
    # xm = read_tsv("species.map", col_names = F)
    # rootnodenumber = 26


    # g2 = g %>% flip(27,31) %>% rotate(21) +

    g2 = g1 %>% rotate(19) %>% rotate(20) %>% rotate(21) %>% rotate(24) + 
        geom_cladelab(node=27, size=1, label="Poaceae", align=TRUE, offset = 0.5, offset.text=0.005) +
        geom_cladelab(node=21, label="Cyperaceae", align=TRUE, offset = 0.5, offset.text=0.005) +
        geom_cladelab(node=26, label="Juncaceae", align=TRUE, offset = 0.5, offset.text=0.005) 


    if (draw == TRUE) {
        pdf(opt$pdf, width = 11, height = 7)
        print(g2)
        dev.off()
        cat("Check: open", opt$pdf, "\n")
    }

    return(g2)
}


# Carex cpdna
step3 <- function(g1,draw=TRUE) {

    g2 = g1 + 
        geom_cladelab(node=34, label="Cyperaceae", align=TRUE, offset = 0.5, offset.text=0.005) +
        geom_cladelab(node=53, label="Juncaceae", align=TRUE, offset = 0.5, offset.text=0.005) 

    if (draw == TRUE) {
        pdf(opt$pdf, width = 11, height = 7)
        print(g2)
        dev.off()
        cat("Check: open", opt$pdf, "\n")
    }

    return(g2)
}

step4 <- function(draw=TRUE) {

    # xm = read_tsv("species.map", col_names = F)
    # xl = read_csv("phylogenies/genes.csv") %>% select(!c("NC_036467","NC_044153","NC_073519"))

    # "chr1_cds-blatx_"
    target_gene = "_cds-blatx_"
    # tl = gsub(target_gene,"", xl$target)
    # tl = gsub("_1","", tl)

    last.target = xl[,ncol(xl)]
    tl = str_extract(last.target[[1]], "blatx_([a-zA-Z0-9]+)_", group=1)
    tl[is.na(tl)] = xl[is.na(tl),ncol(xl)][[1]]


    ngene = nrow(xl)
    nspecies = ncol(xl)
    xl2 = xl
    xl[!is.na(xl)] <- "O"
    xl[is.na(xl)] <- "X"

    dl = data.frame(label=rep(colnames(xl),each=ngene), 
                    category=rep(tl,nspecies), 
                    value=as.numeric(xl=="O"))

    pos.lack.target = c(rep(F,(nspecies-1)*ngene), !grepl(target_gene, last.target))
    dl$value[pos.lack.target] = 0
    dl$value[grepl(",",as.matrix(xl2))] = 2 # no treat for three or more commas 

    dl2 = dl %>% mutate(label2=xm$X2[match(unlist(dl$label), xm$X1)]) %>% select(!label) %>% rename(label=label2)

    p2 <- ggplot(dl2, aes(x=category, y=label, fill=value)) + 
        geom_tile(color="black", show.legend = F) + 
        scale_fill_gradient(low="white", high="grey") + 
        theme(axis.text.x = element_text(angle=90,size=6,hjust = 0.1,vjust=0.2)) + 
        scale_y_discrete(labels = NULL, breaks = NULL) + labs(y = "") +
        xlab(NULL) + ylab(NULL)

    if (draw) {
        pdf(opt$pdf, width = 11, height = 7)
        print(p2)
        dev.off()
        cat("Check: open", opt$pdf, "\n")
    }

    return(p2)
}

# one target
step42 <- function(draw=TRUE) {

    # xm = read_tsv("species.map", col_names = F)
    # xl = read_csv("phylogenies/genes.csv") %>% select(!c("NC_036467","NC_044153","NC_073519"))

    # "chr1_cds-blatx_"
    target_gene = paste0(opt$locus, "_cds-blatx_")

    tl = gsub(target_gene,"", xl$target)
    tl = gsub("_1","", tl)
    ngene = nrow(xl)
    nspecies = ncol(xl)
    xl2 = xl
    xl[!is.na(xl)] <- "O"
    xl[is.na(xl)] <- "X"

    dl = data.frame(label=rep(colnames(xl),each=ngene), 
                    category=rep(tl,nspecies), 
                    value=as.numeric(xl=="O"))

    pos.lack.target = c(rep(F,(nspecies-1)*ngene), !grepl(target_gene, xl2$target))
    dl$value[pos.lack.target] = 0
    dl$value[grepl(",",as.matrix(xl2))] = 2 # no treat for three or more commas 

    dl2 = dl %>% mutate(label2=xm$X2[match(unlist(dl$label), xm$X1)]) %>% select(!label) %>% rename(label=label2)

    p2 <- ggplot(dl2, aes(x=category, y=label, fill=value)) + 
        geom_tile(color="black", show.legend = F) + 
        scale_fill_gradient(low="white", high="grey") + 
        theme(axis.text.x = element_text(angle=90,size=6,hjust = 0.1,vjust=0.2)) + 
        scale_y_discrete(labels = NULL, breaks = NULL) + labs(y = "") +
        xlab(NULL) + ylab(NULL)

    if (draw) {
        pdf(opt$pdf, width = 11, height = 7)
        print(p2)
        dev.off()
        cat("Check: open", opt$pdf, "\n")
    }

    return(p2)
}


# C. pseudochinensis
step5_carex_pseudochinensis <- function(g2,p2,draw=TRUE) {

    # tree = read.newick("phylogenies/s2.treefile", node.label = "support")
    # xm = read_tsv("species.map", col_names = F)
    # rootnodenumber = 26

    pdf(opt$pdf, width = 11, height = 7)
    p2 %>% insert_left(g2, width=1.8) %>% print()
    dev.off()
    cat("Check: open", opt$pdf, "\n")
}

# Carex cpDNA
step5 <- function(g2,p2,draw=TRUE) {

    # tree = read.newick("phylogenies/s2.treefile", node.label = "support")
    # xm = read_tsv("species.map", col_names = F)
    # rootnodenumber = 26

    pdf(opt$pdf, width = 11, height = 7)
    p2 %>% insert_left(g2, width=1.1) %>% print() # CHANGE
    dev.off()
    cat("Check: open", opt$pdf, "\n")
}

# main

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

########################################################################
# main
#
# usage:
#
# /Users/goshng/Dropbox/Documents/Projects/polap/src/run-polap-draw.R --tree phylogenies/s2.treefile --step 1 --pdf t1.pdf
# /Users/goshng/Dropbox/Documents/Projects/polap/src/run-polap-draw.R --map species.map --genes phylogenies/genes.csv --outgroup NC_036467,NC_044153,NC_073519 --step 2 --pdf p2.pdf
# /Users/goshng/Dropbox/Documents/Projects/polap/src/run-polap-draw.R --map species.map --tree phylogenies/s2.treefile --root 26 --step 3 --pdf t2.pdf
# /Users/goshng/Dropbox/Documents/Projects/polap/src/run-polap-draw.R --map species.map --genes phylogenies/genes.csv --outgroup NC_036467,NC_044153,NC_073519 --tree phylogenies/s2.treefile --root 26 --step 4 --pdf tp.pdf

if (opt$step == 1) {
    # /Users/goshng/Dropbox/Documents/Projects/polap/src/run-polap-draw.R --tree phylogenies/s2.treefile --step 1 --pdf t1.pdf
    stopifnot(!is.null(opt$tree))
    step1()
} else if (opt$step == 2) {
    # /Users/goshng/Dropbox/Documents/Projects/polap/src/run-polap-draw.R --map species.map --outgroup NC_036467,NC_044153,NC_073519 --tree phylogenies/s2.treefile --root 26 --step 3 --pdf t2.pdf
    stopifnot(!is.null(opt$tree))
    stopifnot(!is.null(opt$map))
    stopifnot(!is.null(opt$root))
    g = step2()
} else if (opt$step == 3) {
    # /Users/goshng/Dropbox/Documents/Projects/polap/src/run-polap-draw.R --map species.map --tree phylogenies/s2.treefile --root 26 --step 4 --pdf tp.pdf
    stopifnot(!is.null(opt$tree))
    stopifnot(!is.null(opt$map))
    stopifnot(!is.null(opt$root))
    stopifnot(!is.null(opt$node))
    g = step2(FALSE)
    step3(g)
} else if (opt$step == 4) {
    # /Users/goshng/Dropbox/Documents/Projects/polap/src/run-polap-draw.R --map species.map --genes phylogenies/genes.csv --outgroup NC_036467,NC_044153,NC_073519 --step 2 --pdf p2.pdf
    stopifnot(!is.null(opt$map))
    stopifnot(!is.null(opt$genes), !is.null(opt$outgroup))
    p2 = step4()
} else if (opt$step == 5) {
    # /Users/goshng/Dropbox/Documents/Projects/polap/src/run-polap-draw.R --map species.map --genes phylogenies/genes.csv --outgroup NC_036467,NC_044153,NC_073519 --tree phylogenies/s2.treefile --root 26 --step 4 --pdf tp.pdf
    stopifnot(!is.null(opt$tree))
    stopifnot(!is.null(opt$map))
    stopifnot(!is.null(opt$root))
    # stopifnot(!is.null(opt$genes), !is.null(opt$outgroup))
    p2 = step4(FALSE)
    g = step2(FALSE)
    # step4(g)
    g2 = step3(g,FALSE)
    step5(g2, p2)
}



