################################################################################
##### Fig_ones.R
##### Author: Jia Rong Wu, Greg Gloor
##### jwu424 (at) gmail.com
#####
##### DESCRIPTION: Generalized R script in order to generate supporting figures
##### for the paper IQLR. Code for Figure 2 a/b.
#####
##### USAGE: Rscript --vanilla Fig_ones.R
#####
##### LICENSE
##### Copyright (c) 2019 Jia Rong Wu, Greg Gloor
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

################################### FIGURE 4 ###################################
# use reads.1 from Variables.R

conds <- c(rep("A", 10), rep("B", 10))

x.clr.1 <- aldex.clr(reads.1, conds, mc.samples, verbose=FALSE, useMC=FALSE, denom="all")
x.e.1 <- aldex.effect(x.clr.1, useMC=FALSE)
x.t.1  <- aldex.ttest(x.clr.1)
x.all.1 <- data.frame(x.e.1, x.t.1)

x.clr.1 <- aldex.clr(reads.1, conds,  mc.samples, verbose=FALSE, useMC=FALSE, denom="iqlr")
x.e.1 <- aldex.effect(x.clr.1, useMC=FALSE)
x.t.1  <- aldex.ttest(x.clr.1)
x.all.1.iqlr <- data.frame(x.e.1, x.t.1)

x.clr.zero <- aldex.clr(reads.1, conds, mc.samples, verbose=FALSE, denom="zero")
x.e.zero <- aldex.effect(x.clr.zero)
x.t.zero  <- aldex.ttest(x.clr.zero)
x.1.zero <- data.frame(x.e.zero, x.t.zero)

x.clr.lvha <- aldex.clr(reads.1, conds, mc.samples, verbose=FALSE, denom="lvha")
x.e.lvha <- aldex.effect(x.clr.lvha)
x.t.lvha  <- aldex.ttest(x.clr.lvha)
x.1.lvha <- data.frame(x.e.lvha, x.t.lvha)

f1 <- paste(figs.dir, "Fig_ones.pdf",sep="")

plot_ass1 <- function(dataset, main="title"){
	called <- dataset$wi.eBH <= cutoff
	plot(dataset$diff.win, dataset$diff.btw, xlab=xlab, ylab=ylab,
	    col=all.col, pch=all.pch, cex=all.cex, main=main, ylim=c(ymin,ymax))
	points(dataset$diff.win[dataset$rab.all.1 < rare],
	    dataset$diff.btw[dataset$rab.all.1 < rare], col=rare.col, pch=rare.pch,
	    cex=rare.cex)
	points(dataset$diff.win[called], dataset$diff.btw[called], col=called.col,
	    pch=called.pch, cex=called.cex)
	points(dataset[true.set,"diff.win"], dataset[true.set,"diff.btw"],
	    col=true.col, pch=true.pch, cex=true.cex)
	abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
	abline(0,0, col="black")
	abline(0,0, col="white", lwd=thres.lwd, lty=2)

}

pdf(f1, height=9, width=9)
par(fig=c(0,1,0,1), new=TRUE)
ymax=9
par(fig=c(0,0.5,0.5,1), new=TRUE)
    plot_ass1(x.all.1, main="Asymmetric 1")

par(fig=c(0.5,1,0.5,1), new=TRUE)
    plot_ass1(x.all.1.iqlr, main="Asymmetric 1 IQLR")

ymax=8
par(fig=c(0,0.5,0,0.5), new=TRUE)
    plot_ass1(x.1.zero, main="Asymmetric 1 non-zero")

par(fig=c(0.5,1,0,0.5), new=TRUE)
    plot_ass1(x.1.lvha, main="Asymmetric 1 LVHA")


dev.off()

