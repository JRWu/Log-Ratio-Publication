################################################################################ 
##### Analyze.R
##### Author: Jia Rong Wu
##### jwu424 (at) gmail.com 
#####
##### DESCRIPTION: Generalized R script in order to generate supporting figures 
##### for the paper (insert name here) based on the input datasets.
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

##### Check for proper arguments
#args=(commandArgs(trailingOnly=TRUE))
#if (length(args)==0){
#	stop("Expected Usage: Rscript --vanilla Analyze.R filename",call.=FALSE)
#} else if (length(args) == 1){
#	file.name <- as.character(args[1])}

file.name <- "reads_A_500_0"
##### Check if file exists
data.dir <- "../data/"
figs.dir <- "../figures/"

read.file <- paste(data.dir, file.name, sep="")
if (!file.exists(read.file)){
	stop(paste("File: '",file.name,"' does not exist.",sep=""), call.=FALSE)}


##### Load Required Libraries
library(ALDEx2)
library(psych)
library(png)

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
reads[reads==0] <- 0.5
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
thres.line.col="darkgrey"
thres.lwd=1.5





################################### FIGURE 1a ##################################
################################### FIGURE 1a ##################################
################################### FIGURE 1a ##################################
################################### FIGURE 1a ##################################
################################### FIGURE 1a ##################################
reads <- read.table(read.file, header=T, row.names=1, sep="\t", check.names=F)
conds <- c(rep("A", 10), rep("B", 10))
x.clr.o <- aldex.clr(reads, conds, 128, FALSE, FALSE, FALSE)
x.e.o <- aldex.effect(x.clr.o, conds, useMC=TRUE) 
x.t.o  <- aldex.ttest(x.clr.o, conds)
x.all.o <- data.frame(x.e.o, x.t.o)


plot(x.all.o$diff.win, x.all.o$diff.btw, pch=19, col=c(rgb(0,0,0,0.2)), main="Mean Expression Effect Plot",xlab=expression( "Median" ~~ Log[2] ~~ "win-Condition diff" ), ylab=expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" ))
abline(0,0, col=thres.line.col, lty=2, lwd=thres.lwd)



################################### FIGURE 1b ##################################
################################### FIGURE 1b ##################################
################################### FIGURE 1b ##################################
################################### FIGURE 1b ##################################
################################### FIGURE 1b ##################################
dev.new()
hist(x.all.o$diff.btw, breaks=800, xlim=c(-1,1), main="Histogram of Mean Expression", xlab="Mean Expression Per Sample")
abline(v=0, col=thres.line.col, lty=2, lwd=thres.lwd)









################################### FIGURE 2a ##################################
################################### FIGURE 2a ##################################
################################### FIGURE 2a ##################################
################################### FIGURE 2a ##################################
################################### FIGURE 2a ##################################
reads <- read.table(read.file, header=T, row.names=1, sep="\t", check.names=F)
reads <- reads + 0.5

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

dev.new()

png.name <- paste(figs.dir,file.name,"_hist.png",sep="")
png(png.name)
hist(x.all.z$diff.btw, breaks=800, xlim=c(-1,1), main="", xlab="", ylab="")
dev.off()

png.pict <- readPNG(png.name)
aldex.plot(x.all.z)

################################### FIGURE 2b ##################################
################################### FIGURE 2b ##################################
################################### FIGURE 2b ##################################
################################### FIGURE 2b ##################################
################################### FIGURE 2b ##################################
dev.new()
aldex.plot(x.all.o)


################################### FIGURE 3a ##################################
################################### FIGURE 3a ##################################
################################### FIGURE 3a ##################################
################################### FIGURE 3a ##################################
################################### FIGURE 3a ##################################
data(selex)
selex.reads <- selex
selex.conds <- c(rep("NS",7),rep("S",7))

x.clr.o.sel <- aldex.clr(selex.reads, selex.conds, 128, FALSE, FALSE, FALSE)
x.e.o.sel <- aldex.effect(x.clr.o.sel, selex.conds, useMC=TRUE) 
x.t.o.sel  <- aldex.ttest(x.clr.o.sel, selex.conds)
x.all.o.sel <- data.frame(x.e.o.sel, x.t.o.sel)

aldex.plot(x.all.o.sel)
abline(0,0) 



################################### FIGURE 3b ##################################
################################### FIGURE 3b ##################################
################################### FIGURE 3b ##################################
################################### FIGURE 3b ##################################
################################### FIGURE 3b ##################################
##### Does NOT center as well as the zero removal centering
##### Needs to center on the actual data because the invariant set is still influenced by the zero count features

data(selex)
selex.reads <- selex
selex.reads <- selex.reads + 0.5
selex.conds <- c(rep("NS",7),rep("S",7))

selex.reads.clr <- t(apply(selex.reads, 2, function(x){log2(x) - mean(log2(x))}))
selex.reads.var <- apply(selex.reads.clr, 2, function(x){var(x)})
selex.reads.qtl <- quantile(unlist(selex.reads.var))


selex.invariant.set <- which(
	(selex.reads.var < (selex.reads.qtl[upper.bound])) & 
	(selex.reads.var > (selex.reads.qtl[lower.bound]))
)



selex.condition.list <- vector("list", length(unique(selex.conds)))	# list to store conditions
selex.sample.indices <- as.numeric(seq(1, length(selex.conds),1))	# Indicies of samples
selex.feature.indices <- as.numeric(seq(1, nrow(selex.reads), 1))	# Indicies of reads 
selex.neg.indicies <- vector("list", length(unique(selex.conds)))
selex.zero.result <- vector("list", length(unique(selex.conds)))	# list to hold result

for (i in 1:length(unique(selex.conds)))
{
	selex.condition.list[[i]] <- which(selex.conds == unique(selex.conds)[i]) # Condition list
	selex.neg.indicies[[i]] <- selex.invariant.set

}

selex.nr <- nrow( selex.reads )
selex.rn <- rownames( selex.reads )

selex.p <- lapply( selex.reads, function(col) { 
	q <- t( rdirichlet( mc.samples, col)); 
	rownames(q) <- selex.rn; q
})


for (i in 1:length(unique(selex.conds)))
{
	selex.zero.result[[i]] <- lapply( selex.p[selex.condition.list[[i]]], function(m) { 
		apply(log2(m), 2, function(x){mean(x[selex.neg.indicies[[i]]])})
	})
}
selex.set.rev <- unlist(selex.zero.result, recursive=FALSE) # Unlist once to aggregate samples

selex.p.copy <- selex.p
for (i in 1:length(selex.set.rev))
{
	selex.p.copy[[i]] <- as.data.frame(selex.p.copy[[i]])
	selex.p[[i]] <- apply(log2(selex.p.copy[[i]]),1, function(x){ x - (selex.set.rev[[i]])})
	selex.p[[i]] <- t(selex.p[[i]])
}
selex.l2p <- selex.p	# Save the set in order to generate the aldex.clr variable



selex.x <- new("aldex.clr",reads=selex.reads,mc.samples=mc.samples,verbose=verbose,useMC=useMC,analysisData=selex.l2p)

selex.x.tt <- aldex.ttest(selex.x, selex.conds, paired.test=FALSE)
selex.x.effect <- aldex.effect(selex.x, selex.conds, include.sample.summary=include.sample.summary, verbose=verbose)
selex.x.all.z <- data.frame(selex.x.effect, selex.x.tt)

dev.new()
aldex.plot(selex.x.all.z)
abline(0,0)

################################### FIGURE 4a ##################################
################################### FIGURE 4a ##################################
################################### FIGURE 4a ##################################
################################### FIGURE 4a ##################################
################################### FIGURE 4a ##################################
##### Run this set on good data such that it doesn't change











