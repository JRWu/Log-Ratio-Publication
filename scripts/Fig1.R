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
source("Variables.R")

################################### FIGURE 1a ##################################
##### This figure is generated off dataset reads
##### It is the original unmodified base dataset with 40 positive controls
##### This is original ALDEx2 CLR
################################### FIGURE 1a ##################################
##### Read file and ensure it exists
file.name <- "reads"
read.file <- paste(data.dir, file.name, sep="")
if (!file.exists(read.file)){
	stop(paste("File: '",file.name,"' does not exist.",sep=""), call.=FALSE)}

##### Read table and generate conditions
reads <- read.table(read.file, header=T, row.names=1, sep="\t", check.names=F)
x.clr.o <- aldex.clr(reads, conds, mc.samples, zero=FALSE, verbose=FALSE, useMC=FALSE)
conds <- c(rep("A", 10), rep("B", 10))
x.e.o <- aldex.effect(x.clr.o, conds, useMC=TRUE) 
x.t.o  <- aldex.ttest(x.clr.o, conds)
x.all.o <- data.frame(x.e.o, x.t.o)


f1a <- paste(figs.dir, "Fig_1a.png",sep="")
png(f1a)

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
	pushViewport(viewport(x=.2,y=.8,width=.25,height=.25,just=c("left","top")))
	grid.rect()
	par(plt = gridPLT(), new=TRUE)
	hist(x.all.o$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" ), xlab="", ylab="", cex.main=0.8)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)
popViewport(2)

dev.off()

################################### FIGURE 1b ##################################
##### This figure is generated off dataset reads_2
##### It is a dataset with 120 features in condition A that have 0's introduced
##### This is the original ALDEx2 CLR 
################################### FIGURE 1b ##################################
##### Read file and ensure it exists
file.name <- "reads_2"
read.file <- paste(data.dir, file.name, sep="")
if (!file.exists(read.file)){
	stop(paste("File: '",file.name,"' does not exist.",sep=""), call.=FALSE)}

##### Read table and generate conditions
reads <- read.table(read.file, header=T, row.names=1, sep="\t", check.names=F)
conds <- c(rep("A", 10), rep("B", 10))
x.clr.o <- aldex.clr(reads, conds, mc.samples, zero=FALSE, verbose=FALSE, useMC=FALSE)
x.e.o <- aldex.effect(x.clr.o, conds, useMC=TRUE) 
x.t.o  <- aldex.ttest(x.clr.o, conds)
x.all.o <- data.frame(x.e.o, x.t.o)


f1b <- paste(figs.dir, "Fig_1b.png",sep="")
png(f1b)

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
	pushViewport(viewport(x=.2,y=.8,width=.25,height=.25,just=c("left","top")))
	grid.rect()
	par(plt = gridPLT(), new=TRUE)
	hist(x.all.o$diff.btw, breaks=500, xlim=c(-1,1), main=expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" ), xlab="", ylab="", cex.main=0.8)
	abline(v=0, col="black", lwd=thres.lwd)
	abline(v=0, col="white", lwd=thres.lwd, lty=2)
popViewport(2)

dev.off()


