#!/usr/bin/env Rscript

# ./contig-seeds.R contig-annotation-table.txt mt.contig.name-1

suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("readr"))
args = commandArgs(trailingOnly=TRUE)
contig.annotation.table = args[1]
mt.contig.name = args[2]
# contig.annotation.table = "contig-annotation-table.txt"
# mt.contig.name = "mt.contig.name-1"
table.tsv = args[3]

f1 <- function(x) formatC(x, format="fg", big.mark=',')
f2 <- function(x) factor(x, labels=c("", "Yes"))

x = as_tibble(read.table(contig.annotation.table, header = T))
y = as_tibble(read.table(mt.contig.name, sep = '_'))
z = x %>% mutate(Seed = f2(Edge %in% y$V2), L = f1(Length)) %>% select(-Length) %>% relocate(L, .after = Contig) %>% rename(Length=L)
z %>% write.table(sep = '\t', row.names = F, quote = F)
