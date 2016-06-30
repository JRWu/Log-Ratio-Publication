################################################################################ 
##### Fig4.R
##### Author: Jia Rong Wu
##### jwu424 (at) gmail.com 
#####
##### DESCRIPTION: Generalized R script in order to generate supporting figures 
##### for the paper IQLR. Code for Figure 4 a/b.
##### 
##### USAGE: Rscript --vanilla Fig4.R
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

################################### FIGURE 4a ##################################
##### This is a figure generated off dataset count_table2, Tom/Dave's 
##### unpublished dataset.
##### This is the original ALDEx2 transformation.
################################### FIGURE 4a ##################################
file.name <- "count_table2"
read.file <- paste(data.dir, file.name, sep="")
if (!file.exists(read.file)){
	stop(paste("File: '",file.name,"' does not exist.",sep=""), call.=FALSE)}

d <- read.table(read.file, header=T, row.names=1, sep="\t", check.names=F)
reads <- data.frame(t(d["C1_0",6:ncol(d)]),t(d["C2_0",6:ncol(d)]),t(d["C3_0",6:ncol(d)]),t(d["C4_0",6:ncol(d)]),t(d["C5_0",6:ncol(d)]),t(d["WT1_15",6:ncol(d)]),t(d["WT2_15",6:ncol(d)]),t(d["WT3_15",6:ncol(d)]),t(d["WT4_15",6:ncol(d)]),t(d["WT5_15",6:ncol(d)]))
selex.reads <- reads
selex.conds <- c(rep("A",5), rep("B",5))

x.clr.o.sel <- aldex.clr(selex.reads, selex.conds, mc.samples, zero=FALSE, verbose=FALSE, useMC=FALSE)
x.e.o.sel <- aldex.effect(x.clr.o.sel, selex.conds, useMC=TRUE) 
x.t.o.sel  <- aldex.ttest(x.clr.o.sel, selex.conds)
x.all.o.sel <- data.frame(x.e.o.sel, x.t.o.sel)

f4a <- paste(figs.dir, "Fig_4a.png", sep="")
png(f4a)

plot.new()
	pushViewport(viewport())
	called <- x.all.o.sel$we.eBH <= cutoff
	plot(x.all.o.sel$diff.win, x.all.o.sel$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="IQR-Unadjusted Effect Plot")
	points(x.all.o.sel$diff.win[x.all.o.sel$rab.all < rare], x.all.o.sel$diff.btw[x.all.o.sel$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.all.o.sel$diff.win[called], x.all.o.sel$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)
	pushViewport(viewport(x=.65,y=.8,width=.25,height=.25,just=c("left","top")))
	grid.rect()
	par(plt = gridPLT(), new=TRUE)
	hist(x.all.o.sel$diff.btw, breaks=500, xlim=c(-0.5,0.5), main=expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" ), xlab="", ylab="", cex.main=0.8)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)
popViewport(2)

dev.off()

################################### FIGURE 4b ##################################
##### This is a figure generated off dataset count_table2, Tom/Dave's 
##### unpublished dataset.
##### This is the IQLR transformation.
################################### FIGURE 4b ##################################
d <- read.table(read.file, header=T, row.names=1, sep="\t", check.names=F)
reads <- data.frame(t(d["C1_0",6:ncol(d)]),t(d["C2_0",6:ncol(d)]),t(d["C3_0",6:ncol(d)]),t(d["C4_0",6:ncol(d)]),t(d["C5_0",6:ncol(d)]),t(d["WT1_15",6:ncol(d)]),t(d["WT2_15",6:ncol(d)]),t(d["WT3_15",6:ncol(d)]),t(d["WT4_15",6:ncol(d)]),t(d["WT5_15",6:ncol(d)]))
selex.reads <- reads
selex.reads <- selex.reads + 0.5


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


f4b <- paste(figs.dir, "Fig_4b.png", sep="")
png(f4b)


plot.new()
	pushViewport(viewport())
	called <- selex.x.all.z$we.eBH <= cutoff
	plot(selex.x.all.z$diff.win, selex.x.all.z$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="IQR-Adjusted Effect Plot")
	points(selex.x.all.z$diff.win[selex.x.all.z$rab.all < rare], selex.x.all.z$diff.btw[selex.x.all.z$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(selex.x.all.z$diff.win[called], selex.x.all.z$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)
	pushViewport(viewport(x=.65,y=.8,width=.25,height=.25,just=c("left","top")))
	grid.rect()
	par(plt = gridPLT(), new=TRUE)
	hist(selex.x.all.z$diff.btw, breaks=500, xlim=c(-0.5,0.5), main=expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" ), xlab="", ylab="", cex.main=0.8)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)
popViewport(2)

dev.off()
