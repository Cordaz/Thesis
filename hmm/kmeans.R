#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

ncenters <- as.numeric(args[1])
dim <- as.numeric(args[2])
file <- paste(args[3], collapse=NULL)

t <- read.table(file, sep='\t', header=T)

M <- matrix(c(t$SCORE, t$SCORE.1), ncol=2, byrow=F)

km <- kmeans(M, centers=ncenters, iter.max=20, nstart=10)

t$CLUSTER <- km$cluster

library(RColorBrewer)
colors <- brewer.pal(ncenters, "Spectral")

png(paste(as.character(ncenters), "_kmeans.png", collapse=NULL, sep=""), width=7*600, height=7*600, res=600, pointsize=10)
plot(t$SCORE[t$CLUSTER==1], t$SCORE.1[t$CLUSTER==1], col=colors[1], xlim=c(0.5,1), ylim=c(0.5,1), main="KMeans", xlab="Score 0", ylab="Score 1") #yellow
for (i in 2:ncenters) {
	lines(t$SCORE[t$CLUSTER==i], t$SCORE.1[t$CLUSTER==i], col=colors[i], type='p')
}

mu <- matrix(0, ncenters, dim)
sigma <- rep(list(matrix(0, dim, dim)), ncenters)

for(i in 1:ncenters) {
	mu[i,1] <- mean(t$SCORE[t$CLUSTER==i])
	mu[i,2] <- mean(t$SCORE.1[t$CLUSTER==i])
	sd.1 <- sd(t$SCORE[t$CLUSTER==i])
	sd.2 <- sd(t$SCORE.1[t$CLUSTER==i])
	sigma[[i]][1,1] <- sd.1 ^ 2
	sigma[[i]][2,2] <- sd.2 ^ 2
	sigma[[i]][1,2] <- cor(t$SCORE[t$CLUSTER==i], t$SCORE.1[t$CLUSTER==i], method="pearson") * sd.1 * sd.2
	sigma[[i]][2,1] <- cor(t$SCORE.1[t$CLUSTER==i], t$SCORE[t$CLUSTER==i], method="pearson") * sd.1 * sd.2
	#print(sigma[[i]])
}


png(paste(as.character(ncenters), "_emissions_kmeans.png", collapse=NULL, sep=""), width=4*600, height=7*600, res=600, pointsize=10)
par(mfrow=c(ncenters+1,1), mar=c(2,2,1,1))

kd.1 <- density(t$SCORE)
kd.2 <- density(t$SCORE.1)

mean.SCORE <- mean(t$SCORE)
mean.SCORE.1 <- mean(t$SCORE.1)
sd.SCORE <- sd(t$SCORE)
sd.SCORE.1 <- sd(t$SCORE.1)

plot(kd.1$x, kd.1$y * sd.SCORE, col='#ef8a62', type='l', main="Original ditribution", xlab=NULL, ylab=NULL, xlim=c(0.5,1), ylim=c(0, max(max(kd.1$y * sd.SCORE), max(kd.2$y * sd.SCORE.1))))
lines(kd.2$x, kd.2$y * sd.SCORE.1, col='#67a9cf') #blue -> SCORE.1

x <- seq(0.5,1,0.0001)

for (i in 1:ncenters) {
	y.1 <- dnorm(x, mu[i,1], sqrt(sigma[[i]][1,1])) * sqrt(sigma[[i]][1,1])
	#y.1 <- y.1/sum(y.1)
	y.2 <- dnorm(x, mu[i,2], sqrt(sigma[[i]][2,2])) * sqrt(sigma[[i]][1,1])
	#y.2 <- y.2/sum(y.2)
	max.y <- max(y.1, y.2)

	plot(x, y.1, type='l', col='#ef8a62', xlim=c(0.5,1), ylim=c(0,max.y), main=paste("State ", as.character(i-1), collapse=NULL, sep=''))
	lines(x, y.2, type='l', col='#67a9cf')
}

# print(mu)
cat(as.vector(t(mu)), "\n")
for ( i in 1:ncenters) {
	cat(as.vector(t(sigma[[i]])), " ")
}
cat("\n")

dev.off()
