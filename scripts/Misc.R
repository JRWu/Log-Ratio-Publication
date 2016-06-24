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

###########################   FUNCTION_DECLARATIONS   ##########################
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

library(ALDEx2)
library(psych)
library(lattice)
library(gridBase)
library(grid)

figures.dir <- "../figures/"
data.dir <- "../data/"
file.name <- "count_table2"
lower.bound <- 2
upper.bound <- 4
mc.samples <- 128
verbose <- FALSE
useMC <- FALSE
include.sample.summary <- FALSE

f.in <- paste(data.dir, file.name, sep="")


d <- read.table(f.in, header=T, row.names=1, check.names=F, sep="\t")
conds <- c(rep("A",5), rep("B",5))


d.wt.15 <- data.frame(t(d["C1_0",6:ncol(d)]),t(d["C2_0",6:ncol(d)]),t(d["C3_0",6:ncol(d)]),t(d["C4_0",6:ncol(d)]),t(d["C5_0",6:ncol(d)]),t(d["WT1_15",6:ncol(d)]),t(d["WT2_15",6:ncol(d)]),t(d["WT3_15",6:ncol(d)]),t(d["WT4_15",6:ncol(d)]),t(d["WT5_15",6:ncol(d)]))

t.d.clr <- t(apply(d.wt.15, 2, function(x){log2(x) - mean(log2(x))}))
t.d.clr.var <- apply(t.d.clr, 2, function(x){var(x)})
med.t.d.clr <- median(t.d.clr.var)		# Median for the entire set

q <- quantile(unlist(t.d.clr.var))

##### Save the set as a data frame
t.d.clr.var.org <- t.d.clr.var
t.d.clr.var <- as.data.frame(t.d.clr.var)

##### Save the values between Q1 and Q3 
invariant.set <- which(
	(t.d.clr.var[,1] < (q[upper.bound])) & 
	(t.d.clr.var[,1] > (q[lower.bound]))
)

reads <- d.wt.15


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

nr <- nrow( reads )
rn <- rownames( reads )

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






plot.new()
pushViewport(viewport())
called <- x.all.z$we.eBH <= cutoff

plot(x.all.z$diff.win, x.all.z$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="IQR-Adjusted Effect Plot")
abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
points(x.all.z$diff.win[x.all.z$rab.all < rare], x.all.z$diff.btw[x.all.z$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
points(x.all.z$diff.win[called], x.all.z$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
abline(0,0)

pushViewport(viewport(x=.65,y=.8,width=.25,height=.25,just=c("left","top")))
grid.rect()
par(plt = gridPLT(), new=TRUE)
hist(x.all.z$diff.btw, breaks=500, xlim=c(-0.5,0.5), main=expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" ), xlab="", ylab="", cex.main=0.8)
abline(v=0)
popViewport(2)









