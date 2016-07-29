################################################################################
##### Fig1.R
##### Author: Jia Rong Wu
##### jwu424 (at) gmail.com
#####
##### DESCRIPTION: Generalized R script in order to generate supporting figures
##### for the paper IQLR. Code for Figure 1 a/b.
#####
##### USAGE: Rscript --vanilla Fig1.R
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
source("scripts/Variables.R")

################################### FIGURE 1a ##################################
##### This figure is generated off dataset reads
##### It is the original unmodified base dataset with 40 positive controls
##### This is original ALDEx2 denom="all"
################################### FIGURE 1 ##################################
conds <- c(rep("A", 10), rep("B", 10))

##### Read table and generate conditions
x.clr <- aldex.clr(reads, conds, mc.samples, verbose=FALSE, denom="all")
x.e <- aldex.effect(x.clr, conds)
x.t  <- aldex.ttest(x.clr, conds)
x.all <- data.frame(x.e, x.t)

##### Use assymmetric datasets

x.clr.0 <- aldex.clr(reads.0, conds, mc.samples, verbose=FALSE, denom="all")
x.e.0 <- aldex.effect(x.clr.0, conds)
x.t.0  <- aldex.ttest(x.clr.0, conds)
x.all.0 <- data.frame(x.e.0, x.t.0)

x.clr.30 <- aldex.clr(reads.30, conds, mc.samples, verbose=FALSE, denom="all"	)
x.e.30 <- aldex.effect(x.clr.30, conds)
x.t.30  <- aldex.ttest(x.clr.30, conds)
x.all.30 <- data.frame(x.e.30, x.t.30)

f1 <- paste(figs.dir, "Fig_1.pdf",sep="")

pdf(f1, height=6, width=16)
par(fig=c(0,1,0,1), new=TRUE)
par(fig=c(0,0.33,0,1), new=TRUE)
	called <- x.all$wi.eBH <= cutoff
	plot(x.all$diff.win, x.all$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="Symmetric dataset", ylim=c(ymin,ymax))
	points(x.all$diff.win[x.all$rab.all < rare], x.all$diff.btw[x.all$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.all$diff.win[called], x.all$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	points(x.all[true.set,"diff.win"], x.all[true.set,"diff.btw"], col=true.col, pch=true.pch, cex=true.cex)
	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)
par(fig=c(0.04,0.16, 0.4,0.9), new=TRUE)
	hist(x.all$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "btw-Condition diff" ), xlab="", ylab="", cex.main=0.6)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)

par(fig=c(0.33,0.66,0,1), new=TRUE)
	called <- x.all.0$wi.eBH <= cutoff
	plot(x.all.0$diff.win, x.all.0$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="Assymetric 2% dataset", ylim=c(ymin,ymax))
	points(x.all.0$diff.win[x.all.0$rab.all < rare], x.all.0$diff.btw[x.all.0$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.all.0$diff.win[called], x.all.0$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	points(x.all.0[true.set,"diff.win"], x.all.0[true.set,"diff.btw"], col=true.col, pch=true.pch, cex=true.cex)

	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)

par(fig=c(0.37,0.49, 0.4,0.9), new=TRUE)
	hist(x.all.0$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "btw-Condition diff" ), xlab="", ylab="", cex.main=0.6)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)


par(fig=c(0.66,1,0,1), new=TRUE)
	called <- x.all.30$wi.eBH <= cutoff
	plot(x.all.30$diff.win, x.all.30$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="Assymetric 30% dataset", ylim=c(-4,ymax))
	points(x.all.30$diff.win[x.all.30$rab.all < rare], x.all.30$diff.btw[x.all.30$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.all.30$diff.win[called], x.all.30$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	points(x.all.30[true.set,"diff.win"], x.all.30[true.set,"diff.btw"], col=true.col, pch=true.pch, cex=true.cex)

	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)

par(fig=c(0.70,0.84, 0.4,0.9), new=TRUE)
	hist(x.all.30$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "btw-Condition diff" ), xlab="", ylab="", cex.main=0.6)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)

dev.off()


######DESEQ
#condition <- c(rep("A", 10), reap("B", 10))
#cds <- newCountDataSet(reads.30, condition)
# cds = estimateSizeFactors(cds)
# cds <- estimateDispersions(cds, method="per-condition")
#  res <- nbinomTest(cds, "A","B")
#  plot(sqrt(x$fittedDispEsts), res$log2FoldChange)
#  plotDispEsts( cds )
#  plotMA(res)
