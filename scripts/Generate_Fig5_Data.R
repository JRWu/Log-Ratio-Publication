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

data.dir <- "../data/"
figs.dir <- "../figures/"
read.file <- "reads_"

##### Omitting the T-test and only calculating effect sizes reduces mc.samples
mc.samples <- 16
upper.bound <- 4
lower.bound <- 2
verbose <- FALSE
has.BiocParallel <- FALSE
has.parallel <- FALSE
include.sample.summary <- FALSE
useMC <- FALSE

conds <- c(rep("A",10), rep("B",10))

between.mean <- NULL
between.median <- NULL
between.mean.z <- NULL
between.median.z <- NULL
between.mean.o <- NULL
between.median.o <- NULL

##### Iterate through the set of reads with simulated sparsity
for (j in 1:97)
{
	print(j)
	
	##### Zero Transformation
	reads <- read.table(paste(data.dir,read.file,j,sep=""), header=T, row.names=1, sep="\t")
	x.clr.z <- aldex.clr(reads, conds, mc.samples, zero=TRUE, verbose=FALSE)
	x.effect.z <- aldex.effect(x.clr.z, conds, useMC=FALSE, verbose=FALSE) 
	##### Zero Transformation

	##### Original ALDEx2
	reads <- read.table(paste(data.dir,read.file,j,sep=""), header=T, row.names=1, sep="\t")
	x.clr.o <- aldex.clr(reads, conds, mc.samples, zero=FALSE, verbose=FALSE)
	x.effect.o <- aldex.effect(x.clr.o, conds, useMC=FALSE, verbose=FALSE) 
	##### Original ALDEx2

	##### IQR Transformation
	reads <- read.table(paste(data.dir,read.file,j,sep=""), header=T, row.names=1, sep="\t")
	reads <- reads + 0.5
	nr <- nrow( reads )
	rn <- rownames( reads )

	sample.indices <- as.numeric(seq(1, length(conds),1))	# Indicies of samples
	feature.indices <- as.numeric(seq(1, nrow(reads), 1))	# Indicies of reads 
	neg.indicies <- vector("list", length(unique(conds)))
	zero.result <- vector("list", length(unique(conds)))	# list to hold result
	condition.list <- vector("list", length(unique(conds)))	# list to store conditions
	p.by.conditions <- vector("list", length(unique(conds)))

	p.variance <- NULL
	p.variance.avg <- NULL
	p.variance.quantiles <- NULL

	invariant.set <- NULL

	p <- lapply( reads, function(col) { 
		q <- t( rdirichlet( mc.samples, col)); 
		rownames(q) <- rn; q
	})
	
	
    pp <- lapply( p, function(m) {
           apply( log2(m), 2, function(col) { col - mean(col) } )
    })
	
	p.vars <- NULL
	
	##### For every single FEATURE (nr)
	##### Get the feature across all SAMPLES + MC.Instances
	##### 
	
	for (i in 1:nr)
	{
		p.vars <- c(p.vars, var(unlist(lapply(pp,"[",i,))) )
	}
	
	p.vars.quantiles <- quantile(p.vars)

	invariant.set <- which(
		(p.vars < (p.vars.quantiles[upper.bound])) &
		(p.vars > (p.vars.quantiles[lower.bound]))
	)
			
	for (i in 1:length(unique(conds)))
	{
		condition.list[[i]] <- which(conds == unique(conds)[i]) # Condition list
		zero.result[[i]] <- lapply( p[condition.list[[i]]], function(m) { 
			apply(log2(m), 2, function(x){mean(x[invariant.set])})
		})
	}
	set.rev <- unlist(zero.result, recursive=FALSE) # Unlist once to aggregate samples

	p.copy <- p
	for (i in 1:length(set.rev))
	{
		p.copy[[i]] <- as.data.frame(p.copy[[i]])
		p[[i]] <- apply(log2(p.copy[[i]]),1, function(x){ x - (set.rev[[i]])})
		p[[i]] <- t(p[[i]])
	}
	l2p <- p
	
	x <- new("aldex.clr",reads=reads,mc.samples=mc.samples,verbose=verbose,useMC=useMC,analysisData=l2p)

	x.effect <- aldex.effect(x, conds, include.sample.summary=include.sample.summary, verbose=verbose)

	between.mean <- c(between.mean, mean(x.effect$diff.btw))
	between.median <- c(between.median, median(x.effect$diff.btw))	
	
	between.mean.z <- c(between.mean.z, mean(x.effect.z$diff.btw))
	between.median.z <- c(between.median.z, median(x.effect.z$diff.btw))	

	between.mean.o <- c(between.mean.o, mean(x.effect.o$diff.btw))
	between.median.o <- c(between.median.o, median(x.effect.o$diff.btw))	
}

##### Save the 2 sets of data 
f.out.mean <- paste(data.dir,"Instance_Diff_Btw_Mean.txt",sep="")
write.table(between.mean, file=f.out.mean, sep="\t", quote=F)
f.out.median <- paste(data.dir,"Instance_Diff_Btw_Medians.txt",sep="")
write.table(between.median, file=f.out.median, sep="\t", quote=F)

f.out.mean.z <- paste(data.dir,"Zero_Diff_Btw_Mean.txt",sep="")
write.table(between.mean.z, file=f.out.mean.z, sep="\t", quote=F)
f.out.median.z <- paste(data.dir,"Zero_Diff_Btw_Medians.txt",sep="")
write.table(between.median.z, file=f.out.median.z, sep="\t", quote=F)

f.out.mean.o <- paste(data.dir,"Orig_Diff_Btw_Mean.txt",sep="")
write.table(between.mean.o, file=f.out.mean.o, sep="\t", quote=F)
f.out.median.o <- paste(data.dir,"Orig_Diff_Btw_Medians.txt",sep="")
write.table(between.median.o, file=f.out.median.o, sep="\t", quote=F)

