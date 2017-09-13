#!/usr/bin/env Rscript

### parsing arguments and calculating dimensionality
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
	cat("Usage: plot_mix.R <column index (base 1)> <number of mixture components> <input file>\n")
	q()
}
score_col = as.numeric(args[1])
m <- as.numeric(args[2])
file <- paste(args[3], collapse=NULL)

library(mixtools)
library(RColorBrewer)

### reading table
t <- read.table(file, sep='\t', header=T)
xbotlimit <- min(t[score_col])

### performing mixture components identification using EM
mu <- rep(0, m)
sigma <- rep(1, m)
mix <- normalmixEM(t[,score_col], k=m)
for (k in 1:m) {
	mu[k] <- mix$mu[k]
	sigma[k] <- mix$sigma[k]
}

### calculating original distributions
dens <- density(t[, score_col])
kd.x <- dens$x
kd.y <- dens$y
max.kd <- max(kd.y)
t.stddev <- sd(t[, score_col])
x <- seq(xbotlimit,1,0.0001)

### plotting mixture components
png(paste(file, "_mixture.png", collapse=NULL, sep=''), width=4*600, height=5*600, res=600, pointsize=10)
par(mfrow=c(2,1),mar=c(2,2,1,1))

colors <- brewer.pal(m, "Set1")

plot(kd.x, kd.y * t.stddev, col=colors[1], type="l", main="Original ditribution", xlab=NULL, ylab=NULL, xlim=c(xbotlimit,1), ylim=c(0, max.kd * t.stddev))

# genrate values from gaussian and plot it
y <- dnorm(x, mu[1], sigma[1]) * sigma[1]
max.y <- max(y)
for (j in 2:m) {
	y.tmp <- dnorm(x, mu[j], sigma[j]) * sigma[j]
	max.tmp <- max(y)
	if (max.tmp > max.y) {
		max.y <- max.tmp
	}
	y <- cbind(y, y.tmp)
}

plot(x, y[,1], type='l', col=colors[1], xlim=c(xbotlimit,1), ylim=c(0,max.y), main="Mixture component")
for (j in 2:m) {
	lines(x, y[,j], type='l', col=colors[j])
}

#Corner_text <- function(text, location="topright"){
#	legend(location,legend=text, bty ="n", pch=NA)
#}
#mix_infos = "Comp.1: mu = "
#mix_infos <- paste(mix_infos, as.character(mu[1]), sep="", collapse=NULL)
#mix_infos <- paste(mix_infos, "\tsigma^2 = ", sep="", collapse=NULL)
#mix_infos <- paste(mix_infos, as.character(sigma[1]), sep="", collapse=NULL)
#mix_infos = paste(mix_infos, "\nComp.2: mu = ", sep="", collapse=NULL)
#mix_infos <- paste(mix_infos, as.character(mu[2]), sep="", collapse=NULL)
#mix_infos <- paste(mix_infos, "\tsigma^2 = ", sep="", collapse=NULL)
#mix_infos <- paste(mix_infos, as.character(sigma[2]), sep="", collapse=NULL)
#Corner_text(text=mix_infos,location= "bottomright")

dev.off()
