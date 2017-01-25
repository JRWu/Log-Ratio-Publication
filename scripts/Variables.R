################################################################################
##### Variables.R
##### Author: Jia Rong Wu
##### jwu424 (at) gmail.com
#####
##### DESCRIPTION: R script that loads variables and plotting variables common
##### to each figures script for breveity.
#####
##### USAGE: At the beginning of each script that requires these variables,
##### append the following line: source("Variables.R")
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

##### Setup correct directories
data.dir <- "data/"
figs.dir <- "figures/"

##### Load Required Libraries
library(ALDEx2)
#library(psych)
#library(lattice)
#library(gridBase)
#library(grid)

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

# this is the dataset
# can easily remake the raw data by introducing 0 or 1 count variables in
# one or other subset at runtime
##### variables 47-56 5X
##### variables 57-66 4X
##### variables 67-76 3X
##### variables 77-86 2X
##### variable 108 is FP effect >1, variables 1679, 3671 fp effect < -1:  just over threshold by chance
file.name <- "d.n0.txt" # this is the original reads file with 0 count variables removed
read.file <- paste(data.dir, file.name, sep="")

if (!file.exists(read.file)){
	stop(paste("File: '",file.name,"' does not exist.",sep=""), call.=FALSE)}

reads.all <- read.table(read.file, header=T, row.names=1, sep="\t", check.names=F)

# reduce to first 1000 lines
reads <- reads.all[1:1000,]

for(i in seq(47,86, by=2)){
	r <- as.numeric(c(reads[i,11:20], reads[i,1:10]))
	reads[i,] <- r
}



true.set <- rownames(reads[47:86,])


# reads with 2% of genes with large assymmetry as 0 values
reads.0 <- reads
reads.0[980:1000,1:10] <- 0

# reads with 6% assymmetry as 0 values
reads.30 <- reads
reads.30[940:1000,1:10] <- 0


# reads with 2% assymmetry as 1 values
reads.1 <- reads
reads.1[980:1000,1:10] <- 1

# reads with 2% assymmetry as 10 values
reads.10 <- reads
reads.10[980:1000,1:10] <- 10

# reads with 2% assymmetry as 0 values, orthogonal to real difference
reads.orth <- reads
reads.orth[940:1000,5:15] <- 0

# reads with 2% assymmetry as 0 values, 45 degrees to real difference
reads.diag <- reads
reads.diag[940:1000,3:13] <- 0


upper.bound <- 4
lower.bound <- 2
verbose <- FALSE
useMC <- FALSE
include.sample.summary <- FALSE
mc.samples <- 16

##### Plotting Variable Declarations
xlab=NULL
ylab=NULL
xlim=NULL
ylim=NULL
all.col=rgb(0,0,0,0.2)
all.pch=19
all.cex=0.4
called.col="red"
called.pch=20
called.cex=0.6
true.col="blue"
true.pch=1
true.cex=0.7
thres.line.col="darkgrey"
thres.lwd=1.5
test="wilcox"
cutoff=0.1
rare.col="black"
rare=0
rare.pch=20
rare.cex=0.2
xlab <- expression( "Median" ~~ Log[2] ~~ "dispersion" )
ylab <- expression( "Median" ~~ Log[2] ~~ "difference" )
ymin <- -3
ymax <- 12
