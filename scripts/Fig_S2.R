################################################################################ 
##### Fig_S2.R
##### Author: Jia Rong Wu
##### jwu424 (at) gmail.com 
#####
##### DESCRIPTION: Generalized R script in order to generate supporting figures 
##### for the paper IQLR. Generates figures that depict what an effect plot is 
##### relative to a boxplot.
##### 
##### USAGE: Rscript --vanilla Fig_S2.R
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

################################### FIGURE S2a #################################
##### This figure is generated off dataset reads
##### It is the original unmodified base dataset with 40 positive controls
##### This is original ALDEx2 CLR
################################### FIGURE S2a #################################
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

s2a <- paste(figs.dir, "Fig_S2a.png", sep="")
png(s2a)

called <- x.all.o$we.eBH <= cutoff
plot(x.all.o$diff.win, x.all.o$diff.btw, xlab=xlab, ylab=ylab, col=all.col, pch=all.pch, cex=all.cex, main="ALDEx2 Effect Plot", ylim=c(ymin,ymax))
points(x.all.o$diff.win[x.all.o$rab.all < rare], x.all.o$diff.btw[x.all.o$rab.all < rare], col=rare.col, pch=rare.pch, cex=rare.cex)
points(x.all.o$diff.win[called], x.all.o$diff.btw[called], col=called.col, pch=called.pch, cex=called.cex)
abline(0,1, col=thres.line.col, lty=2, lwd=thres.lwd)
abline(0,-1, col=thres.line.col, lty=2, lwd=thres.lwd)
abline(0,0, col="black")
abline(0,0, col="white", lwd=thres.lwd, lty=2)

dev.off()


################################### FIGURE S2b #################################
##### This figure depicts a boxplot with error bars.
################################### FIGURE S2b #################################

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
	if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
	stop("vectors must be same length")
	arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ..., lwd=2)
}

y <- rnorm(100, mean=0.6, sd=0.4)
y.2 <- rnorm(100, mean=1.1, sd=1)
y <- c(y,y.2)
y <- matrix(y,100,2)
y.means <- apply(y,2,mean)
y.sd <- apply(y,2,sd)



s2b <- paste(figs.dir, "Fig_S2b.png", sep="")
png(s2b)

barx <- barplot(y.means, names.arg=c("A","B"),ylim=c(0,1.5), col=c(rgb(0,0,0,0.4)), axis.lty=1, xlab="Conditions", ylab="Value (arbitrary units)", main="Barplot", cex.main=2)
error.bar(barx,y.means, 1.96*y.sd/10)

dev.off()
