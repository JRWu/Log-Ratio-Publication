################################################################################ 
##### Fig3.R
##### Author: Jia Rong Wu
##### jwu424 (at) gmail.com 
#####
##### DESCRIPTION: Generalized R script in order to generate supporting figures 
##### for the paper IQLR. Code for Figure 3 a/b.
##### 
##### USAGE: Rscript --vanilla Fig3.R 
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

##### Load Required Libraries
source("Variables.R")

################################### FIGURE 3a ##################################
##### This figure is generated off dataset reads_2_min
##### It is a dataset where 120 features have their values substituted for 1's
##### This is the zero-removal transformation
################################### FIGURE 3a ##################################
##### Read file and ensure it exists
file.name <- "reads_2_min"
read.file <- paste(data.dir, file.name, sep="")
if (!file.exists(read.file)){
	stop(paste("File: '",file.name,"' does not exist.",sep=""), call.=FALSE)}

reads <- read.table(read.file, header=T, row.names=1, sep="\t", check.names=F)
conds <- c(rep("A", 10), rep("B", 10))

x.clr.o <- aldex.clr(reads, conds, mc.samples , zero=TRUE, FALSE, FALSE)
x.e.o <- aldex.effect(x.clr.o, conds, useMC=TRUE) 
x.t.o  <- aldex.ttest(x.clr.o, conds)
x.all.o <- data.frame(x.e.o, x.t.o)

f3a <- paste(figs.dir, "Fig_3a.png", sep="")
png(f3a)

plot.new()
	pushViewport(viewport())
	called <- x.all.o$we.eBH <= cutoff
	plot(x.all.o$diff.win, x.all.o$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="IQR-Unadjusted Effect Plot", ylim=c(ymin,ymax))
	points(x.all.o$diff.win[x.all.o$rab.all < rare], x.all.o$diff.btw[x.all.o$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.all.o$diff.win[called], x.all.o$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)
	pushViewport(viewport(x=.65,y=.8,width=.25,height=.25,just=c("left","top")))
	grid.rect()
	par(plt = gridPLT(), new=TRUE)
	hist(x.all.o$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" ), xlab="", ylab="", cex.main=0.8)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)
popViewport(2)

dev.off()

################################### FIGURE 3b ##################################
##### This figure is generated off dataset reads_2_min
##### It is a dataset where 120 features have their values substituted for 1's
##### This is the IQLR transformation
################################### FIGURE 3b ##################################
##### Read file and ensure it exists
file.name <- "reads_2_min"
read.file <- paste(data.dir, file.name, sep="")
if (!file.exists(read.file)){
	stop(paste("File: '",file.name,"' does not exist.",sep=""), call.=FALSE)}

reads <- read.table(read.file, header=T, row.names=1, sep="\t", check.names=F)
conds <- c(rep("A", 10), rep("B", 10))

#####  remove all rows with reads less than the minimum set by minsum 
minsum <- 0

##### remove any row in which the sum of the row is 0
z <- as.numeric(apply(reads, 1, sum))
reads <- as.data.frame( reads[(which(z > minsum)),]  )

##### Adjust all reads with prior of 0.5
reads <- reads + 0.5

nr <- nrow( reads )
rn <- rownames( reads )

##### Generate the CLR transformation of the DATA and get variance of the CLR
reads.clr <- t(apply(reads, 2, function(x){log2(x) - mean(log2(x))}))
reads.var <- apply(reads.clr, 2, function(x){var(x)})
reads.qtl <- quantile(unlist(reads.var))

##### Get the indicies of the "invariant set" features
invariant.set <- which(
	(reads.var < (reads.qtl[upper.bound])) & 
	(reads.var > (reads.qtl[lower.bound]))
)

##### Setup general lists in order to store the intermediate values
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


##### Generate the Geometric Mean on the invariant set(s)
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


##### Generate the aldex.clr object on the IQLR transformed reads
x <- new("aldex.clr",reads=reads,mc.samples=mc.samples,verbose=verbose,useMC=useMC,analysisData=l2p)

x.tt <- aldex.ttest(x, conds, paired.test=FALSE)
x.effect <- aldex.effect(x, conds, include.sample.summary=include.sample.summary, verbose=verbose)
x.all.z <- data.frame(x.effect, x.tt)

f3b <- paste(figs.dir, "Fig_3b.png", sep="")
png(f3b)

plot.new()
	pushViewport(viewport())
	called <- x.all.z$we.eBH <= cutoff
	plot(x.all.z$diff.win, x.all.z$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="IQR-Adjusted Effect Plot", ylim=c(ymin,ymax))
	points(x.all.z$diff.win[x.all.z$rab.all < rare], x.all.z$diff.btw[x.all.z$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.all.z$diff.win[called], x.all.z$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)
	
	pushViewport(viewport(x=.65,y=.8,width=.25,height=.25,just=c("left","top")))
	grid.rect()
	par(plt = gridPLT(), new=TRUE)
	hist(x.all.z$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" ), xlab="", ylab="", cex.main=0.8)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)
popViewport(2)

dev.off()
