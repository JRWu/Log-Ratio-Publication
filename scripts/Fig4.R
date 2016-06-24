################################################################################ 
##### Fig4.R
##### Author: Jia Rong Wu
##### jwu424 (at) gmail.com 
#####
##### DESCRIPTION: Generalized R script in order to generate supporting figures 
##### for the paper (insert name here). Code for Figure 4 a/b.
##### 
##### USAGE: Rscript --vanilla Analyze.R filename
##### 
##### LICENSE 
##### Copyright (c) 2016 Jia Rong Wu
##### 
##### Permission is hereby granted, free of charge, to any person obtaining a 
##### copy of this software and associated documentation files (the "Software"), 
##### to deal in the Software without restriction, including without limitation 
##### the rights to use, copy, modify, merge, publish, distribute, sublicense, 
##### and/or sell copies of the Software, and to permit persons to whom the 
##### Software is furnished to do so, subject to the following conditions:
##### 
##### The above copyright notice and this permission notice shall be included in 
##### all copies or substantial portions of the Software.
##### 
##### THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
##### IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
##### FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
##### THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
##### LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
##### FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
##### DEALINGS IN THE SOFTWARE.
################################################################################ 

file.name <- "reads"
data.dir <- "../data/"
figs.dir <- "../figures/"

read.file <- paste(data.dir, file.name, sep="")
if (!file.exists(read.file)){
	stop(paste("File: '",file.name,"' does not exist.",sep=""), call.=FALSE)}


##### Load Required Libraries
library(ALDEx2)
library(psych)
library(lattice)
library(gridBase)
library(grid)

##### Function Declarations
rdirichlet <- function (n, alpha)
{
  if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) can not be negative in rtriang")
  if(is.vector(alpha)) alpha <- t(alpha)
  l <- dim(alpha)[2]
  x <- matrix(rgamma(l * n, t(alpha)), ncol = l, byrow=TRUE)  # Gere le recycling
  return(x / rowSums(x))
}

#################### Read Tables and Variable Declarations  ####################
reads <- read.table(read.file, header=T, row.names=1, sep="\t", check.names=F)
conds <- c(rep("A", 10), rep("B", 10))

upper.bound <- 4
lower.bound <- 2

verbose <- FALSE
has.BiocParallel <- FALSE
has.parallel <- FALSE
include.sample.summary <- FALSE
useMC <- FALSE
zero <- TRUE
mc.samples <- 128
as.numeric(as.integer(mc.samples))
minsum <- 0
prior <- 0.5
nr <- nrow( reads )
rn <- rownames( reads )

xlab=NULL
ylab=NULL
xlim=NULL
ylim=NULL
all.col=rgb(0,0,0,0.2)
all.pch=19
all.cex=0.4
called.col="red"
called.pch=20
called.cex=0.6
thres.line.col="darkgrey"
thres.lwd=1.5
test="welch"
cutoff=0.1
rare.col="black"
rare=0
rare.pch=20
rare.cex=0.2
xlab <- expression( "Median" ~~ Log[2] ~~ "win-Condition diff" )
ylab <- expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" )



################################### FIGURE 4a ##################################
################################### FIGURE 4a ##################################
################################### FIGURE 4a ##################################
################################### FIGURE 4a ##################################
################################### FIGURE 4a ##################################
reads <- read.table(read.file, header=T, row.names=1, sep="\t", check.names=F)
conds <- c(rep("A", 10), rep("B", 10))
x.clr.o <- aldex.clr(reads, conds, 128, FALSE, FALSE, FALSE)
x.e.o <- aldex.effect(x.clr.o, conds, useMC=TRUE) 
x.t.o  <- aldex.ttest(x.clr.o, conds)
x.all.o <- data.frame(x.e.o, x.t.o)





f2a <- paste(figs.dir,"Fig_4a_IQR_Unadjusted.png",sep="")
png(f2a)

plot.new()
pushViewport(viewport())
called <- x.all.o$we.eBH <= cutoff

plot(x.all.o$diff.win, x.all.o$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="IQR-Unadjusted Effect Plot")
points(x.all.o$diff.win[x.all.o$rab.all < rare], x.all.o$diff.btw[x.all.o$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
points(x.all.o$diff.win[called], x.all.o$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
abline(0,0)

pushViewport(viewport(x=.65,y=.8,width=.25,height=.25,just=c("left","top")))
grid.rect()
par(plt = gridPLT(), new=TRUE)
hist(x.all.o$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" ), xlab="", ylab="", cex.main=0.8)
abline(v=0)
popViewport(2)

dev.off()


################################### FIGURE 4b ##################################
################################### FIGURE 4b ##################################
################################### FIGURE 4b ##################################
################################### FIGURE 4b ##################################
################################### FIGURE 4b ##################################
reads <- read.table(read.file, header=T, row.names=1, sep="\t", check.names=F)

#####  remove all rows with reads less than the minimum set by minsum 
minsum <- 0

##### remove any row in which the sum of the row is 0
z <- as.numeric(apply(reads, 1, sum))
reads <- as.data.frame( reads[(which(z > minsum)),]  )
reads <- reads + 0.5

nr <- nrow( reads )
rn <- rownames( reads )


reads.clr <- t(apply(reads, 2, function(x){log2(x) - mean(log2(x))}))
reads.var <- apply(reads.clr, 2, function(x){var(x)})
reads.qtl <- quantile(unlist(reads.var))

invariant.set <- which(
	(reads.var < (reads.qtl[upper.bound])) & 
	(reads.var > (reads.qtl[lower.bound]))
)

condition.list <- vector("list", length(unique(conds)))	# list to store conditions
sample.indices <- as.numeric(seq(1, length(conds),1))	# Indicies of samples
feature.indices <- as.numeric(seq(1, nrow(reads), 1))	# Indicies of reads 
neg.indicies <- vector("list", length(unique(conds)))
zero.result <- vector("list", length(unique(conds)))	# list to hold result

for (i in 1:length(unique(conds)))
{
	condition.list[[i]] <- which(conds == unique(conds)[i]) # Condition list
	neg.indicies[[i]] <- invariant.set

}

p <- lapply( reads, function(col) { 
	q <- t( rdirichlet( mc.samples, col)); 
	rownames(q) <- rn; q
})


for (i in 1:length(unique(conds)))
{
	zero.result[[i]] <- lapply( p[condition.list[[i]]], function(m) { 
		apply(log2(m), 2, function(x){mean(x[neg.indicies[[i]]])})
	})
}
set.rev <- unlist(zero.result, recursive=FALSE) # Unlist once to aggregate samples

p.copy <- p
for (i in 1:length(set.rev))
{
	p.copy[[i]] <- as.data.frame(p.copy[[i]])
	p[[i]] <- apply(log2(p.copy[[i]]),1, function(x){ x - (set.rev[[i]])})
	p[[i]] <- t(p[[i]])
}
l2p <- p	# Save the set in order to generate the aldex.clr variable



x <- new("aldex.clr",reads=reads,mc.samples=mc.samples,verbose=verbose,useMC=useMC,analysisData=l2p)

x.tt <- aldex.ttest(x, conds, paired.test=FALSE)
x.effect <- aldex.effect(x, conds, include.sample.summary=include.sample.summary, verbose=verbose)
x.all.z <- data.frame(x.effect, x.tt)




f2b <- paste(figs.dir,"Fig_4b_IQR_Adjusted.png",sep="")
png(f2b)
plot.new()
pushViewport(viewport())
called <- x.all.z$we.eBH <= cutoff

plot(x.all.z$diff.win, x.all.z$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="IQR-Adjusted Effect Plot")
points(x.all.z$diff.win[x.all.z$rab.all < rare], x.all.z$diff.btw[x.all.z$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
points(x.all.z$diff.win[called], x.all.z$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
abline(0,0)

pushViewport(viewport(x=.65,y=.8,width=.25,height=.25,just=c("left","top")))
grid.rect()
par(plt = gridPLT(), new=TRUE)
hist(x.all.z$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" ), xlab="", ylab="", cex.main=0.8)
abline(v=0)
popViewport(2)

dev.off()









































