################################################################################ 
##### Fig1.R
##### Author: Jia Rong Wu
##### jwu424 (at) gmail.com 
#####
##### DESCRIPTION: Generalized R script in order to generate supporting figures 
##### for the paper (insert name here). Code for Figure 1 a/b.
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

file.name <- "reads_A_500_0"
data.dir <- "../data/"
figs.dir <- "../figures/"

read.file <- paste(data.dir, file.name, sep="")
if (!file.exists(read.file)){
	stop(paste("File: '",file.name,"' does not exist.",sep=""), call.=FALSE)}


##### Load Required Libraries
library(ALDEx2)
library(psych)
library(png)

##### Function Declarations
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

#################### Read Tables and Variable Declarations  ####################
reads <- read.table(read.file, header=T, row.names=1, sep="\t", check.names=F)
reads[reads==0] <- 0.5
conds <- c(rep("A", 10), rep("B", 10))

upper.bound <- 4
lower.bound <- 2

verbose <- FALSE
has.BiocParallel <- FALSE
has.parallel <- FALSE
include.sample.summary <- FALSE
useMC <- FALSE
zero <- TRUE
mc.samples <- 128
as.numeric(as.integer(mc.samples))
minsum <- 0
prior <- 0.5
nr <- nrow( reads )
rn <- rownames( reads )
thres.line.col="darkgrey"
thres.lwd=1.5


################################### FIGURE 1a ##################################
################################### FIGURE 1a ##################################
################################### FIGURE 1a ##################################
################################### FIGURE 1a ##################################
################################### FIGURE 1a ##################################
reads <- read.table(read.file, header=T, row.names=1, sep="\t", check.names=F)
conds <- c(rep("A", 10), rep("B", 10))
x.clr.o <- aldex.clr(reads, conds, 128, FALSE, FALSE, FALSE)
x.e.o <- aldex.effect(x.clr.o, conds, useMC=TRUE) 
x.t.o  <- aldex.ttest(x.clr.o, conds)
x.all.o <- data.frame(x.e.o, x.t.o)


f1a <- paste(figs.dir,"Fig_1a_EffectPlot.png", sep="")

png(f1a)
plot(x.all.o$diff.win, x.all.o$diff.btw, pch=19, col=c(rgb(0,0,0,0.2)), main="Mean Expression Effect Plot",xlab=expression( "Median" ~~ Log[2] ~~ "win-Condition diff" ), ylab=expression( "Median" ~~ Log[2] ~~ "btw-Condition diff" ))
abline(0,0, col=thres.line.col, lty=2, lwd=thres.lwd)
dev.off()


################################### FIGURE 1b ##################################
################################### FIGURE 1b ##################################
################################### FIGURE 1b ##################################
################################### FIGURE 1b ##################################
################################### FIGURE 1b ##################################

f1b <- paste(figs.dir,"Fig_1b_EffectHistogram.png", sep="")

png(f1b)
hist(x.all.o$diff.btw, breaks=800, xlim=c(-1,1), main="Histogram of Mean Expression", xlab="Mean Expression Per Sample")
abline(v=0, col=thres.line.col, lty=2, lwd=thres.lwd)
dev.off()









