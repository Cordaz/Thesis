#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
file = paste(args, collapse=NULL)
#file
t <- read.table(file, sep="\t", header=TRUE)
t_sorted <- t[order(t$diff_log2, decreasing=TRUE), ]

out_f = paste(file, 'sorted', sep='.', collapse=NULL)

write.table(t_sorted, out_f, sep='\t', append=FALSE, quote=FALSE, eol='\n', row.names=FALSE, col.names=names(t_sorted))
