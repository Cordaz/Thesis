#!/usr/bin/env Rscript

### parsing arguments and calculating dimensionality
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
	cat("Usage: plot_mix.R <comma-separate list of column index (base 1)> <number of mixture components> <input file>\n")
	cat("Number of states is determined as N_mix_comp ^ N_features\n")
	q()
}
scores_list = strsplit(paste(args[1], collapse=NULL), ",")
score_col = c()
for(i in 1:length(scores_list)) {
	score_col <- c(score_col, as.numeric(scores_list[[i]]))
}
nfeatures <- length(score_col)
m <- as.numeric(args[2])
file <- paste(args[3], collapse=NULL)
nstates <- m^nfeatures

library(mixtools)
library(RColorBrewer)

### combine function
# creates all combination of a [nxm] matrix to a [m^n, n] matrix, with one element per row
# | a b c d |    | a e i |
# | e f g h | -> | a e j |
# | i j k l \    |  ...  |
#                | a h l |
#                | b e i |
#                |  ...  |
#                | d h l |

combine <- function(M, row, C, n, m, start) {
	# check if going down out of M
	if (row > dim(M)[1]) {
		return (C)
	}
	slot <- m ^ (n - 1)

	for (x in 1:m) {
		for (i in 1:slot) {
			C[start + (x-1)*slot + (i-1), row] <- M[row,x]
		}
		C <- combine(M, row+1, C, n-1, m, start + (x-1)*slot)
	}
	return (C)
}


### reading table
t <- read.table(file, sep='\t', header=T)

### performing mixture components identification using EM
means <- matrix(0, nrow=nfeatures, ncol=m)
vars <- matrix(1, nrow=nfeatures, ncol=m)
j <- 1
for (i in score_col) {
	mix <- normalmixEM(t[,i], k=m)
	for (k in 1:m) {
		means[j,k] <- mix$mu[k]
		vars[j,k] <- mix$sigma[k]
	}
	j <- j + 1
}

### create all combinations of state using the combine function
mu <- matrix(0, nstates, nfeatures)
sigma <- matrix(0, nstates, nfeatures)

mu <- combine(means, 1, mu, nfeatures, m, 1)
sigma <- combine(vars, 1, sigma, nfeatures, m, 1)

### calculating original distributions
t.stddev <- rep(0, nfeatures)
dens <- density(t[, score_col[1]])
kd.x <- dens$x
kd.y <- dens$y
t.stddev[1] <- sd(t[, score_col[1]])
max.kd <- max(kd.y) * t.stddev[1]
for (i in 2:nfeatures) {
	dens <- density(t[, score_col[i]])
	kd.x <- cbind(kd.x, dens$x)
	kd.y <- cbind(kd.y, dens$y)
	t.stddev[i] <- sd(t[, score_col[i]])
	if (max(dens$y) * t.stddev[i] > max.kd) {
		max.kd <- max(dens$y) * t.stddev[i]
	}
}
x <- seq(0.5,1,0.0001)

### plotting mixture components
png(paste(as.character(nstates), "_mixture.png", collapse=NULL, sep=''), width=4*600, height=5*600, res=600, pointsize=10)
par(mfrow=c(nfeatures+1,1),mar=c(2,2,1,1))

colors <- brewer.pal(nfeatures, "Set1")

plot(kd.x[,1], kd.y[,1] * t.stddev[1], col=colors[1], type="l", main="Original ditribution", xlab=NULL, ylab=NULL, xlim=c(0.5,1), ylim=c(0, max.kd))
for (i in 2:nfeatures) {
	lines(kd.x[,i], kd.y[,i] * t.stddev[i], col=colors[i])
}

# genrate values from gaussian and plot it
for (i in 1:nfeatures) {
	y <- dnorm(x, means[i,1], vars[i,1]) * vars[i,1]
	max.y <- max(y)
	for (j in 2:m) {
		y.tmp <- dnorm(x, means[i,j], vars[i,j]) * vars[i,j]
		max.tmp <- max(y)
		if (max.tmp > max.y) {
			max.y <- max.tmp
		}
		y <- cbind(y, y.tmp)
	}

	plot(x, y[,1], type='l', col=colors[i], xlim=c(0.5,1), ylim=c(0,max.y), main=paste("Feature ", as.character(i), collapse=NULL, sep=''))
	for (j in 2:m) {
		lines(x, y[,j], type='l', col=colors[i])
	}
}

dev.off()

### assign states
STATES <- rep(0, dim(t)[1])
PROB_STATES <- matrix(0, nrow=dim(t)[1], ncol=nstates)

for (i in 1:dim(t)[1]) { # for each sample
	maxp <- 0
	for (s in 1:nstates) { # for each state
		p <- 1
		for (j in 1:length(score_col)) { # for each feature
			p <- p * dnorm(t[i,score_col[j]], mu[s,j], sigma[s,j]) * sigma[s,j]
		}
		PROB_STATES[i,s] <- p
		if (p > maxp) {
			maxp <- p
			STATES[i] <- s - 1
		}
	}
}

t <- cbind(t, STATES)
names <- paste("P_STATE_0", collapse=NULL, sep='')
for (i in 2:nstates) {
	names <- cbind(names, paste("P_STATE_", as.character(i-1), collapse=NULL, sep=''))
}
colnames(PROB_STATES) <- names
t <- cbind(t, PROB_STATES)


# cannot do in more than 2 dim
if (nfeatures == 2) {
	colors <- brewer.pal(nstates, "Set1")
	png(paste(as.character(nstates), "_states_mix.png", collapse=NULL, sep=""), width=7*600, height=7*600, res=600, pointsize=10)
	plot(t[t$STATES==0,score_col[1]], t[t$STATES==0, score_col[2]], col=colors[1], xlim=c(0.5,1), ylim=c(0.5,1), main="Assigned states", xlab="Score 0", ylab="Score 1")
	for (i in 2:nstates) {
		lines(t[t$STATES==i-1, score_col[1]], t[t$STATES==i-1, score_col[2]], col=colors[i], type='p')
	}

	dev.off()
}

counts <- rep(0, nstates)
for (i in 1:dim(t)[1]) {
	for (j in 1:nstates) {
		if (t$STATES[i] == j-1) {
			counts[j] <- counts[j] + 1
		}
	}
}

### creating zscore of means
t.means <- rep(0, nfeatures)
for (i in 1:length(score_col)) {
	t.means[i] <- mean(t[, score_col[i]])
}
mu.zscore <- matrix(0, nstates, nfeatures)
mu.zscore <- (mu - t.means) / ( t.stddev / sqrt(counts) )
max.muzscore <- max(mu.zscore)
min.muzscore <- min(mu.zscore)

### creating states graph

png(paste(as.character(nstates), "_emissions_mix.png", collapse=NULL, sep=""), width=4*600, height=7*600, res=600, pointsize=10)
par(mfrow=c(nstates+1,1), mar=c(2,2,1,1))

dev.base <- dev.cur()

png(paste(as.character(nstates), "_emissions_mix_zscore.png", collapse=NULL, sep=""), width=4*600, height=7*600, res=600, pointsize=10)
par(mfrow=c(nstates,1), mar=c(2,2,1,1))

dev.zscore <- dev.cur()

dev.set(dev.base)
plot(kd.x[,1], kd.y[,1] * t.stddev[1], col=colors[1], type="l", main="Original ditribution", xlab=NULL, ylab=NULL, xlim=c(0.5,1), ylim=c(0, max.kd))
for (i in 2:nfeatures) {
	lines(kd.x[,i], kd.y[,i] * t.stddev[i], col=colors[i])
}

for (i in 1:nstates) {
	y <- dnorm(x, mu[i,1], sigma[i,1]) * sigma[i,1]
	max.y <- max(y)
	for (j in 2:nfeatures) {
		y.tmp <- dnorm(x, mu[i,j], sigma[i,j]) * sigma[i,j]
		max.tmp <- max(y)
		if (max.tmp > max.y) {
			max.y <- max.tmp
		}
		y <- cbind(y, y.tmp)
	}

	dev.set(dev.base)
	plot(x, y[,1], type='l', col=colors[1], xlim=c(0.5,1), ylim=c(0,max.y), main=paste("State ", as.character(i-1), collapse=NULL, sep=''))
	for (j in 2:nfeatures) {
		lines(x, y[,j], type='l', col=colors[j])
	}

	dev.set(dev.zscore)
	plot(rep(mu.zscore[i,1], 2), c(0,1), type='l', col=colors[1], xlim=c(min.muzscore, max.muzscore), ylim=c(0,1), main=paste("State ", as.character(i-1), collapse=NULL, sep=''))
	for (j in 2:nfeatures) {
		lines(rep(mu.zscore[i,j], 2), c(0,1), col=colors[j])
	}
}

dev.off(dev.base)
dev.off(dev.zscore)


### print output
# create names
names <- "Component 1"
for (i in 2:nfeatures) {
	names <- cbind(names, paste("Component ", as.character(i), collapse=NULL, sep=''))
}
colnames(mu) <- names
colnames(mu.zscore) <- names
colnames(sigma) <- names

names <- "State 0"
for (i in 2:nstates) {
	names <- cbind(names, paste("State ", as.character(i-1), collapse=NULL, sep=''))
}
rownames(mu) <- names
rownames(mu.zscore) <- names
rownames(sigma) <- names
names(counts) <- names
outfile <- paste(file, paste(as.character(nstates), "states.stats", collapse=NULL,sep=''), collapse=NULL, sep='')
cat(as.character(nstates), file=outfile, append=F)
cat(" states: stats\n>Means\n\t\t", file=outfile, append=T)
write.table(mu, file=outfile, append=T, sep='\t', eol='\n', quote=F)
cat("\n>Means (Z-score)\n\t\t", file=outfile, append=T)
write.table(mu.zscore, file=outfile, append=T, sep='\t', eol='\n', quote=F)
cat("\n>Variances\n\t\t", file=outfile, append=T)
write.table(sigma, file=outfile, append=T, sep='\t', eol='\n', quote=F)
cat("\n>Counts\n", file=outfile, append=T)
write.table(counts, file=outfile, append=T, sep='\t', eol='\n', quote=F, col.names=F)

write.table(t, paste(file, paste(as.character(nstates), "states.tables", collapse=NULL,sep=''), collapse=NULL, sep=''), sep="\t", append=FALSE, quote=FALSE, eol="\n", row.names=FALSE, col.names=names(t))
