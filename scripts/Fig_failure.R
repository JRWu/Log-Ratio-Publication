################################################################################
##### Fig5.R
##### Author: Jia Rong Wu
##### jwu424 (at) gmail.com
#####
##### DESCRIPTION: Generalized R script in order to generate supporting figures
##### for the paper IQLR. Code for Figure 5.
#####
##### USAGE: Rscript --vanilla Fig5.R
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

################################### FIGURE 5 ###################################
##### This figure is generated off data from the script Generate_Fig5_Data.R
##### This dataset demonstrates how each transformation is affected by
##### systematic variation.
################################### FIGURE 5 ###################################
data.dir <- "data/"
##### Read IQR Transformation Data
instance.median <- unlist(read.table(file=paste(data.dir,"Instance_Diff_Btw_Medians.txt",sep=""), header=T, row.names=1, sep="\t"))
instance.mean <- unlist(read.table(file=paste(data.dir,"Instance_Diff_Btw_Mean.txt",sep=""), header=T, row.names=1, sep="\t"))

##### Read Zero Transformation Data
zero.median <- unlist(read.table(file=paste(data.dir,"Zero_Diff_Btw_Medians.txt",sep=""), header=T, row.names=1, sep="\t"))
zero.mean <- unlist(read.table(file=paste(data.dir,"Zero_Diff_Btw_Mean.txt",sep=""), header=T, row.names=1, sep="\t"))

##### Read LVHA Transformation Data
lvha.median <- unlist(read.table(file=paste(data.dir,"LVHA_Btw_Medians.txt",sep=""), header=T, row.names=1, sep="\t"))
lvha.mean <- unlist(read.table(file=paste(data.dir,"LVHA_Btw_Mean.txt",sep=""), header=T, row.names=1, sep="\t"))

##### Read Median Transformation Data
med.median <- unlist(read.table(file=paste(data.dir,"Median_Diff_Btw_Medians.txt",sep=""), header=T, row.names=1, sep="\t"))
med.mean <- unlist(read.table(file=paste(data.dir,"Median_Diff_Btw_Mean.txt",sep=""), header=T, row.names=1, sep="\t"))

##### Read Original Transformation Data
original.median <- unlist(read.table(file=paste(data.dir,"Orig_Diff_Btw_Medians.txt",sep=""), header=T, row.names=1, sep="\t"))
original.mean <- unlist(read.table(file=paste(data.dir,"Orig_Diff_Btw_Mean.txt",sep=""), header=T, row.names=1, sep="\t"))

# des.median <- unlist(read.table(file=paste(data.dir,"DES_Btw_Medians.txt",sep=""), # header=T, row.names=1, sep="\t"))

##### Plot Work
ymax <- max(zero.median, instance.median, original.median)
ymin <- min(zero.median, instance.median, original.median)

ylb <- "Median Between Condition Diff"
xlb <- "Percent Asymmetric Sparsity"
man <- paste(ylb, "  vs. ", xlb, sep="")

darkness <- 0.5
pch <- c(15,16,17,19,5,19)
col <- c( rgb(0,0,1,darkness),rgb(1,0,0,darkness), rgb(0,1,1,darkness),rgb(1,0.3,0,1),rgb(0,0,0,1))
leg <- c("CLR","NZLR","IQLR","MED","LVHA")

fig.name <- paste(figs.dir,"Fig_failure.pdf",sep="")
pdf(fig.name, height=5, width=8)

##### Generate an Empty Plot
plot(NULL,xlim=c(0,50),ylim=c(ymin,ymax), xlab=xlb, ylab=ylb) #, main=man)

##### Plot the points per transformation
points(zero.median[1:50], pch=pch[2], col=col[2], type="b")
points(original.median[1:50], pch=pch[1], col=col[1], type="b")
points(med.median[1:50], pch=pch[4], col=col[4], type="b")
points(lvha.median[1:50], pch=pch[5], col=col[5], type="b",cex=1.2)
points(instance.median[1:50], pch=pch[3], col=col[3], type="b")
#points(des.median[1:50], pch=pch[6], col="black", type="b",cex=1.2)

abline(0,0, col=c(rgb(0,0,0,0.75)))
abline(v=25, col=c(rgb(0,0,0,0.2)))
abline(v=50, col=c(rgb(0,0,0,0.2)))

##### Position the legend at the top right of the plot
legend(0,7,legend=leg, pch=pch, col=col)

dev.off()
