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

x.clr.o.sel <- aldex.clr(selex.reads, mc.samples, verbose=FALSE, useMC=FALSE)
x.e.o.sel <- aldex.effect(x.clr.o.sel, selex.conds, useMC=TRUE)
x.t.o.sel  <- aldex.ttest(x.clr.o.sel, selex.conds)
x.all.o.sel <- data.frame(x.e.o.sel, x.t.o.sel)

x.clr.osel.iqlr <- aldex.clr(selex.reads, mc.samples, verbose=FALSE, useMC=FALSE, iqlr=TRUE)
x.e.osel.iqlr <- aldex.effect(x.clr.osel.iqlr, selex.conds, useMC=TRUE)
x.t.osel.iqlr  <- aldex.ttest(x.clr.osel.iqlr, selex.conds)
x.all.osel.iqlr <- data.frame(x.e.osel.iqlr, x.t.osel.iqlr)

f4 <- paste(figs.dir, "Fig_4.pdf", sep="")
ymax=2
ymin= -0.5
pdf(f4, height=6, width=11)
par(fig=c(0,1,0,1), new=TRUE)
par(fig=c(0,0.5,0,1), new=TRUE)
	called <- x.all.o.sel$wi.eBH <= cutoff
	plot(x.all.o.sel$diff.win, x.all.o.sel$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="SELEX", ylim=c(ymin,ymax))
	points(x.all.o.sel$diff.win[x.all.o.sel$rab.all.iqlr < rare], x.all.o.sel$diff.btw[x.all.o.sel$rab.all.iqlr < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.all.o.sel$diff.win[called], x.all.o.sel$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	points(x.all.o.sel$diff.win[called], x.all.o.sel$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)
#par(fig=c(0.05,0.25, 0.4,0.9), new=TRUE)
#	hist(x.all.o.sel$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "btw-Condition diff" ), xlab="", ylab="", cex.main=0.6)
#	abline(v=0, col="black", lwd=thres.lwd)
#	abline(v=0, col="white", lwd=thres.lwd, lty=2)

################################### FIGURE 1b ##################################
##### This figure is generated off dataset reads_2
##### It is a dataset with 120 features in condition A that have 0's introduced
##### This is the original ALDEx2 CLR
################################### FIGURE 1b ##################################

par(fig=c(0.5,1,0,1), new=TRUE)
	called <- x.all.osel.iqlr$wi.eBH <= cutoff
	plot(x.all.osel.iqlr$diff.win, x.all.osel.iqlr$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="SELEX IQLR-adjusted", ylim=c(ymin,ymax))
	points(x.all.osel.iqlr$diff.win[x.all.osel.iqlr$rab.all < rare], x.all.osel.iqlr$diff.btw[x.all.osel.iqlr$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
	points(x.all.osel.iqlr$diff.win[called], x.all.osel.iqlr$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
	#points(x.all.osel.iqlr[invariant.set,"diff.win"], x.all.osel.iqlr[invariant.set,"diff.btw"], col=rgb(0,0,1,0.1), pch=called.pch, cex=called.cex)

	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)

#par(fig=c(0.55,0.75, 0.4,0.9), new=TRUE)
#	hist(x.all.osel.iqlr$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "btw-Condition diff" ), xlab="", ylab="", cex.main=0.6)
#	abline(v=0, col="black", lwd=thres.lwd)
#	abline(v=0, col="white", lwd=thres.lwd, lty=2)

dev.off()
