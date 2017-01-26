################################################################################
##### Fig2.R
##### Author: Jia Rong Wu
##### jwu424 (at) gmail.com
#####
##### DESCRIPTION: Generalized R script in order to generate supporting figures
##### for the paper IQLR. Code for Figure 2 a/b.
#####
##### USAGE: Rscript --vanilla Fig2.R
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

################################### FIGURE 2 ###################################
##### This figure is generated off dataset reads_2
##### It is a dataset where 120 features have their values substituted for 0's
##### This is the IQLR transformation that FIXES the problem
##### requires new ALDEx2
################################### FIGURE 2 ###################################
# use reads.1 and reads.10 from Variables.R

x.clr.1 <- aldex.clr(reads.1, mc.samples, verbose=FALSE, useMC=FALSE, iqlr=FALSE)
x.e.1 <- aldex.effect(x.clr.1, conds, useMC=TRUE)
x.t.1  <- aldex.ttest(x.clr.1, conds)
x.all.1 <- data.frame(x.e.1, x.t.1)

x.clr.10 <- aldex.clr(reads.10, mc.samples, verbose=FALSE, useMC=FALSE, iqlr=FALSE)
x.e.10 <- aldex.effect(x.clr.10, conds, useMC=TRUE)
x.t.10  <- aldex.ttest(x.clr.10, conds)
x.all.10 <- data.frame(x.e.10, x.t.10)


x.clr.1 <- aldex.clr(reads.1, mc.samples, verbose=FALSE, useMC=FALSE, iqlr=TRUE)
x.e.1 <- aldex.effect(x.clr.1, conds, useMC=TRUE)
x.t.1  <- aldex.ttest(x.clr.1, conds)
x.all.1.iqlr <- data.frame(x.e.1, x.t.1)

x.clr.10 <- aldex.clr(reads.10, mc.samples, verbose=FALSE, useMC=FALSE, iqlr=TRUE)
x.e.10 <- aldex.effect(x.clr.10, conds, useMC=TRUE)
x.t.10  <- aldex.ttest(x.clr.10, conds)
x.all.10.iqlr <- data.frame(x.e.10, x.t.10)


f1 <- paste(figs.dir, "Fig_3.pdf",sep="")

pdf(f1, height=9, width=9)
par(fig=c(0,1,0,1), new=TRUE)
ymax=9
par(fig=c(0,0.5,0.5,1), new=TRUE)
	called <- x.all.1$wi.eBH <= cutoff
	plot(x.all.1$diff.win, x.all.1$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="Aymmetric 1", ylim=c(ymin,ymax))
	points(x.all.1$diff.win[x.all.1$rab.all.1 < rare], x.all.1$diff.btw[x.all.1$rab.all.1 < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.all.1$diff.win[called], x.all.1$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	points(x.all.1$diff.win[called], x.all.1$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)

par(fig=c(0.5,1,0.5,1), new=TRUE)
	called <- x.all.1.iqlr$wi.eBH <= cutoff
	plot(x.all.1.iqlr$diff.win, x.all.1.iqlr$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="Symmetric 1 IQLR-adjusted", ylim=c(ymin,ymax))
	points(x.all.1.iqlr$diff.win[x.all.1.iqlr$rab.all.1.iqlr < rare], x.all.1.iqlr$diff.btw[x.all.1.iqlr$rab.all.1.iqlr < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.all.1.iqlr$diff.win[called], x.all.1.iqlr$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	points(x.all.1.iqlr$diff.win[called], x.all.1.iqlr$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)

################################### FIGURE 1b ##################################
##### This figure is generated off dataset reads_2
##### It is a dataset with 120 features in condition A that have 0's introduced
##### This is the original ALDEx2 CLR
################################### FIGURE 1b ##################################
ymax=7
par(fig=c(0,0.5,0,0.5), new=TRUE)
	called <- x.all.10$wi.eBH <= cutoff
	plot(x.all.10$diff.win, x.all.10$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="Aymmetric 10", ylim=c(ymin,ymax))
	points(x.all.10$diff.win[x.all.10$rab.all.10 < rare], x.all.10$diff.btw[x.all.10$rab.all.10 < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.all.10$diff.win[called], x.all.10$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	points(x.all.10$diff.win[called], x.all.10$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)

par(fig=c(0.5,1,0,0.5), new=TRUE)
	called <- x.all.10.iqlr$wi.eBH <= cutoff
	plot(x.all.10.iqlr$diff.win, x.all.10.iqlr$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="Aymmetric 10 IQLR-adjusted", ylim=c(ymin,ymax))
	points(x.all.10.iqlr$diff.win[x.all.10.iqlr$rab.all.10.iqlr < rare], x.all.10.iqlr$diff.btw[x.all.10.iqlr$rab.all.10.iqlr < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.all.10.iqlr$diff.win[called], x.all.10.iqlr$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	points(x.all.10.iqlr$diff.win[called], x.all.10.iqlr$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)


dev.off()

