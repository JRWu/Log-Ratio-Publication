################################################################################
##### Fig1.R
##### Author: Jia Rong Wu, Greg Gloor
##### jwu424 (at) gmail.com
#####
##### DESCRIPTION: Generalized R script in order to generate supporting figures
##### for the paper IQLR. Code for Figure 2.
#####
##### USAGE: Rscript --vanilla Fig1.R
#####
##### LICENSE
##### Copyright (c) 2016, 2019 Jia Rong Wu, GBG
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
mc.samples=128
################################### FIGURE 1a ##################################
##### This figure is generated off dataset reads
##### It is the original unmodified base dataset with 40 positive controls
##### This is ALDEx2 with denom="iqlr | zero | custom"
################################### FIGURE 1 ##################################

##### Read table and generate conditions
conds <- c(rep("A", 10), rep("B", 10))

##### Use assymmetric datasets

x.clr.iqlr <- aldex.clr(reads.0, conds, mc.samples, verbose=FALSE, denom="iqlr")
x.e.iqlr <- aldex.effect(x.clr.iqlr)
x.t.iqlr  <- aldex.ttest(x.clr.iqlr)
x.iqlr.iqlr <- data.frame(x.e.iqlr, x.t.iqlr)


x.clr.zero <- aldex.clr(reads.0, conds, mc.samples, verbose=FALSE, denom="zero")
x.e.zero <- aldex.effect(x.clr.zero)
x.t.zero  <- aldex.ttest(x.clr.zero)
x.zero.zero <- data.frame(x.e.zero, x.t.zero)
zero.feature.subset <- aldex.set.mode(reads.0, conds, denom="zero")


x.clr.user <- aldex.clr(reads.0, conds, mc.samples, verbose=FALSE, denom="lvha")
x.e.user <- aldex.effect(x.clr.user)
x.t.user  <- aldex.ttest(x.clr.user)
x.user.user <- data.frame(x.e.user, x.t.user)
denom=c(100:200)

x.clr.med <- aldex.clr(reads.0, conds, mc.samples, verbose=FALSE, denom="median")
x.e.med <- aldex.effect(x.clr.med)
x.t.med  <- aldex.ttest(x.clr.med)
x.med.med <- data.frame(x.e.med, x.t.med)

f2 <- paste(figs.dir, "Fig_2.pdf",sep="")

pdf(f2, height=9, width=9)
par(fig=c(0,1,0,1), new=TRUE)

par(fig=c(0,0.5,0.5,1), new=TRUE)
	called <- x.iqlr.iqlr$wi.eBH <= cutoff
	plot(x.iqlr.iqlr$diff.win, x.iqlr.iqlr$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="IQLR set", ylim=c(ymin,ymax))
	points(x.iqlr.iqlr$diff.win[x.iqlr.iqlr$rab.all < rare], x.iqlr.iqlr$diff.btw[x.iqlr.iqlr$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.iqlr.iqlr$diff.win[called], x.iqlr.iqlr$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	points(x.iqlr.iqlr[true.set,"diff.win"], x.iqlr.iqlr[true.set,"diff.btw"], col=true.col, pch=true.pch, cex=true.cex)
	points(x.iqlr.iqlr[x.clr.iqlr@denom,"diff.win"], x.iqlr.iqlr[x.clr.iqlr@denom,"diff.btw"], col=rgb(0.8,0.7,0.5,0.5), pch=true.pch, cex=true.cex)

	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)

par(fig=c(0.05,0.3, 0.7,0.95), new=TRUE)
	hist(x.iqlr.iqlr$diff.btw, breaks=500, xlim=c(-1,1), main=NULL, xlab="", ylab="", cex.main=0.5, cex.lab=0.5)
	title("Difference", line=0.5, cex.main=0.6)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)

par(fig=c(0.5,1,0.5,1), new=TRUE)
	called <- x.zero.zero$wi.eBH <= cutoff
	plot(x.zero.zero$diff.win, x.zero.zero$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="non-zero set", ylim=c(ymin,ymax))
	points(x.zero.zero$diff.win[x.zero.zero$rab.all < rare], x.zero.zero$diff.btw[x.zero.zero$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.zero.zero$diff.win[called], x.zero.zero$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	points(x.zero.zero[true.set,"diff.win"], x.zero.zero[true.set,"diff.btw"], col=true.col, pch=true.pch, cex=true.cex)
	points(x.zero.zero[x.clr.zero@denom,"diff.win"], x.zero.zero[x.clr.zero@denom,"diff.btw"], col=rgb(0.8,0.7,0.5,0.5), pch=true.pch, cex=true.cex)

	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)

par(fig=c(0.55,0.8, 0.7,0.95), new=TRUE)
	hist(x.zero.zero$diff.btw, breaks=500, xlim=c(-1,1), main=NULL, xlab="", ylab="", cex.main=0.5)
	title("Difference", line=0.5, cex.main=0.6)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)

par(fig=c(0,0.5,0,0.5), new=TRUE)
	called <- x.user.user$wi.eBH <= cutoff
	plot(x.user.user$diff.win, x.user.user$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="LVHA set", ylim=c(ymin,ymax))
	points(x.user.user$diff.win[x.user.user$rab.all < rare], x.user.user$diff.btw[x.user.user$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.user.user$diff.win[called], x.user.user$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	points(x.user.user[true.set,"diff.win"], x.user.user[true.set,"diff.btw"], col=true.col, pch=true.pch, cex=true.cex)
	points(x.user.user[x.clr.user@denom,"diff.win"], x.user.user[x.clr.user@denom,"diff.btw"], col=rgb(0.8,0.7,0.5,0.5), pch=true.pch, cex=true.cex)

	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)

par(fig=c(0.05,0.3,0.2,0.45), new=TRUE)
	hist(x.user.user$diff.btw, breaks=500, xlim=c(-1,1), main=NULL, xlab="", ylab="", cex.main=0.5)
	title("Difference", line=0.5, cex.main=0.6)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)

par(fig=c(0.5,1,0,0.5), new=TRUE)
	called <- x.med.med$wi.eBH <= cutoff
	plot(x.med.med$diff.win, x.med.med$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="median", ylim=c(ymin,ymax))
	points(x.med.med$diff.win[x.med.med$rab.all < rare], x.med.med$diff.btw[x.med.med$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.med.med$diff.win[called], x.med.med$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	points(x.med.med[true.set,"diff.win"], x.med.med[true.set,"diff.btw"], col=true.col, pch=true.pch, cex=true.cex)

	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)

par(fig=c(0.55,0.8,0.2,0.45), new=TRUE)
	hist(x.med.med$diff.btw, breaks=500, xlim=c(-1,1), main=NULL, xlab="", ylab="", cex.main=0.5)
	title("Difference", line=0.5, cex.main=0.6)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)

dev.off()

#######DESEQ
#condition <- c(rep("A", 10), reap("B", 10))
#cds <- newCountDataSet(reads.30, condition)
# cds = estimateSizeFactors(cds)
# cds <- estimateDispersions(cds, method="per-condition")
#  res <- nbinomTest(cds, "A","B")
#  plot(sqrt(res$fittedDispEsts), res$log2FoldChange)
#  plotDispEsts( cds )
#  plotMA(res)
