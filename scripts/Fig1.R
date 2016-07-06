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
##### This is original ALDEx2 CLR
################################### FIGURE 1a ##################################

##### Read table and generate conditions
x.clr <- aldex.clr(reads, mc.samples, verbose=FALSE, useMC=FALSE)
conds <- c(rep("A", 10), rep("B", 10))
x.e <- aldex.effect(x.clr, conds, useMC=TRUE)
x.t  <- aldex.ttest(x.clr, conds)
x.all <- data.frame(x.e, x.t)

##### make assymmetric dataset
reads.0 <- reads
reads.0[5800:5859,1:10] <- 0

x.clr.0 <- aldex.clr(reads.0, mc.samples, verbose=FALSE, useMC=FALSE)
x.e.0 <- aldex.effect(x.clr.0, conds, useMC=TRUE)
x.t.0  <- aldex.ttest(x.clr.0, conds)
x.all.0 <- data.frame(x.e.0, x.t.0)


f1 <- paste(figs.dir, "Fig_1.pdf",sep="")

pdf(f1, height=6, width=11)
par(fig=c(0,1,0,1), new=TRUE)
par(fig=c(0,0.5,0,1), new=TRUE)
	called <- x.all$we.eBH <= cutoff
	plot(x.all$diff.win, x.all$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="IQR-Unadjusted Effect Plot", ylim=c(ymin,ymax))
	points(x.all$diff.win[x.all$rab.all < rare], x.all$diff.btw[x.all$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.all$diff.win[called], x.all$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	points(x.all$diff.win[called], x.all$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)
par(fig=c(0.05,0.25, 0.4,0.9), new=TRUE)
	hist(x.all$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" ), xlab="", ylab="", cex.main=0.6)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)

################################### FIGURE 1b ##################################
##### This figure is generated off dataset reads_2
##### It is a dataset with 120 features in condition A that have 0's introduced
##### This is the original ALDEx2 CLR
################################### FIGURE 1b ##################################

par(fig=c(0.5,1,0,1), new=TRUE)
	called <- x.all.0$we.eBH <= cutoff
	plot(x.all.0$diff.win, x.all.0$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="IQR-Unadjusted Effect Plot", ylim=c(ymin,ymax))
	points(x.all.0$diff.win[x.all.0$rab.all < rare], x.all.0$diff.btw[x.all.0$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.all.0$diff.win[called], x.all.0$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)

par(fig=c(0.55,0.75, 0.4,0.9), new=TRUE)
	hist(x.all.0$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" ), xlab="", ylab="", cex.main=0.6)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)
#popViewport(2)

dev.off()


