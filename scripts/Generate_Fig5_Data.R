################################################################################
##### Generate_Fig5_Data.R
##### Author: Jia Rong Wu
##### jwu424 (at) gmail.com
#####
##### DESCRIPTION: Generalized R script in order to generate supporting figures
##### for the paper (insert name here). Generates data for comparison purposes.
##### Runs slowly. Expect it to take ~2-3 hours on a mobile class i7 processor.
#####
##### USAGE: Rscript --vanilla Generate_Fig5_Data.R
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
library(ALDEx2)

data.dir <- "data/"
figs.dir <- "figures/"
read.file <- "sparse/reads_"

##### Omitting the T-test and only calculating effect sizes reduces mc.samples
mc.samples <- 64
upper.bound <- 4
lower.bound <- 2
verbose <- FALSE
has.BiocParallel <- TRUE
has.parallel <- TRUE
include.sample.summary <- FALSE
useMC <- FALSE

conds <- c(rep("A",10), rep("B",10))

between.mean.i <- NULL
between.median.i <- NULL
between.mean.z <- NULL
between.median.z <- NULL
between.mean.o <- NULL
between.median.o <- NULL
between.mean.m <- NULL
between.median.m <- NULL
between.mean.l <- NULL
between.median.l <- NULL

##### Iterate through the set of reads with simulated sparsity
for (j in 1:24)
{
	print(j)
	reads <- read.table(paste(data.dir,read.file,j,sep=""), header=T, row.names=1, sep="\t")

	##### LVHA Transformation
	x.clr.l <- aldex.clr(reads, conds, mc.samples, denom="lvha", verbose=FALSE, useMC=FALSE)
	x.effect.l <- aldex.effect(x.clr.i, conds, useMC=FALSE, verbose=FALSE)

#	##### Zero Transformation
#	x.clr.z <- aldex.clr(reads, conds, mc.samples, denom="zero", verbose=FALSE,useMC=FALSE)
#	x.effect.z <- aldex.effect(x.clr.z, conds, useMC=FALSE, verbose=FALSE)
#	##### Zero Transformation
#
#	##### Original ALDEx2
#	x.clr.o <- aldex.clr(reads, conds, mc.samples, denom="all", verbose=FALSE, useMC=FALSE)
#	x.effect.o <- aldex.effect(x.clr.o, conds, useMC=FALSE, verbose=FALSE)
#	##### Original ALDEx2
#
#	##### IQR Transformation
#	x.clr.i <- aldex.clr(reads, conds, mc.samples, denom="iqlr", verbose=FALSE, useMC=FALSE)
#	x.effect.i<- aldex.effect(x.clr.i, conds, useMC=FALSE, verbose=FALSE)
#
#	##### median Transformation
#	x.clr.m <- aldex.clr(reads, conds, mc.samples, denom="median", verbose=FALSE, useMC=FALSE)
#	x.effect.m<- aldex.effect(x.clr.m, conds, useMC=FALSE, verbose=FALSE)
#
#
#	between.mean.i <- c(between.mean.i, mean(x.effect.i$diff.btw))
#	between.median.i <- c(between.median.i, median(x.effect.i$diff.btw))
#
	between.mean.l <- c(between.mean.l, mean(x.effect.l$diff.btw))
	between.median.l <- c(between.median.l, median(x.effect.l$diff.btw))
#
#	between.mean.z <- c(between.mean.z, mean(x.effect.z$diff.btw))
#	between.median.z <- c(between.median.z, median(x.effect.z$diff.btw))
#
#	between.mean.o <- c(between.mean.o, mean(x.effect.o$diff.btw))
#	between.median.o <- c(between.median.o, median(x.effect.o$diff.btw))
#
#	between.mean.m <- c(between.mean.m, mean(x.effect.m$diff.btw))
#	between.median.m <- c(between.median.m, median(x.effect.m$diff.btw))
}

##### Save the 2 sets of data
f.out.mean <- paste(data.dir,"LVHA_Btw_Mean.txt",sep="")
write.table(between.mean.l, file=f.out.mean, sep="\t", quote=F)
f.out.median <- paste(data.dir,"LVHA_Btw_Medians.txt",sep="")
write.table(between.median.l, file=f.out.median, sep="\t", quote=F)

#f.out.mean <- paste(data.dir,"Instance_Diff_Btw_Mean.txt",sep="")
#write.table(between.mean.i, file=f.out.mean, sep="\t", quote=F)
#f.out.median <- paste(data.dir,"Instance_Diff_Btw_Medians.txt",sep="")
#write.table(between.median.i, file=f.out.median, sep="\t", quote=F)
#
#f.out.mean.z <- paste(data.dir,"Zero_Diff_Btw_Mean.txt",sep="")
#write.table(between.mean.z, file=f.out.mean.z, sep="\t", quote=F)
#f.out.median.z <- paste(data.dir,"Zero_Diff_Btw_Medians.txt",sep="")
#write.table(between.median.z, file=f.out.median.z, sep="\t", quote=F)
#
#f.out.mean.o <- paste(data.dir,"Orig_Diff_Btw_Mean.txt",sep="")
#write.table(between.mean.o, file=f.out.mean.o, sep="\t", quote=F)
#f.out.median.o <- paste(data.dir,"Orig_Diff_Btw_Medians.txt",sep="")
#write.table(between.median.o, file=f.out.median.o, sep="\t", quote=F)
#
#f.mut.mean.m <- paste(data.dir,"Median_Diff_Btw_Mean.txt",sep="")
#write.table(between.mean.m, file=f.mut.mean.m, sep="\t", quote=F)
#f.mut.median.m <- paste(data.dir,"Median_Diff_Btw_Medians.txt",sep="")
#write.table(between.median.m, file=f.mut.median.m, sep="\t", quote=F)
#
