#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE);
file <- args[1]

# Load packages
# Check if installed
list.of.packages <- c("seqLogo")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) {
	source("http://bioconductor.org/biocLite.R")
	biocLite("seqLogo")
}
library(seqLogo)

png(paste(file, ".png", collapse=NULL, sep=''), width=8*600, height=3*600, res=600, pointsize=10)

t <- read.table(file, sep='\t', header=F)
M <- as.matrix(t(t))
colnames(M) <- c('A','C','G','T')
rownames(M) <- NULL
df <- as.data.frame(M)

#define function that divides the frequency by the row sum i.e. proportions
proportion <- function(x){
   rs <- sum(x)
   return(x / rs)
}

#create position weight matrix
pwm <- apply(df, 1, proportion)
pwm <- makePWM(pwm)
seqLogo(pwm)
dev.off()
