################################################################################ 
##### Fig_S1.R
##### Author: Jia Rong Wu
##### jwu424 (at) gmail.com 
#####
##### DESCRIPTION: Generalized R script in order to generate supporting figures 
##### for the paper IQLR. Generates figures that depict how the 
##### Geometric mean is influenced by sparse data. 
##### 
##### USAGE: Rscript --vanilla Fig_S1.R
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
library(psych)

################################### FIGURE S1a #################################
##### Generate a dataset of 150 random variates. 10% of these variates are 0's.
##### Used to illustrate the difference in geometric means by computing with the
##### entire set of 150 variates, or only the variates that fall within the IQR 
##### of the dataset. 
################################### FIGURE S1a #################################

figs.dir <- "../figures/"
gm.name <- paste(figs.dir, "Fig_S1a.png",sep="")
png(gm.name)			# Save file in figures directory

##### Generate randomly drawn data 135 variates, mean of 20, std. dev of 3
x <- rnorm(n=135, m= 20, sd=3)
prior <- c(rep(0.5,15))
combined <- c(x,prior)				# Simulate 10% zero features in 150 numbers

list_histo <- hist(combined, col=c(rgb(0,0,0,0.15)), main="Geometric Means", xlab="Count", freq=TRUE,cex.main=2.0, breaks=15)

len <-max(list_histo$counts) + 1	# Describes how many values are in histogram (used for letter position and arrow positioning)
clip(0,100,0,100)					# Restricts the ablines to graph limits
arrow.len <- 2.5					# For length of the arrows
upper.bound <- 4					# Upper bound of Quantile 3
lower.bound <- 2					# Lower bound of Quantile 2

q.combined <- quantile(combined)	# Draw quantiles

##### NOTE: quantiles are formatted as such:
#####      0%      25%      50%      75%     100% 
##### 0.50000 17.26034 19.38153 21.07459 28.48562 
##### Therefore q.combined[4] refers to the upper 75%
##### and 		q.combined[2] refers to the lower 25%

##### Find the invariant set that fall within quantiles 2 and 3 
combined.invariant.set <- which(
	(combined < q.combined[upper.bound]) &
	(combined > q.combined[lower.bound])
)

##### Compute geometric means of both adjusted and non-adjusted variables
org.gmean <- geometric.mean(combined)
iqr.gmean <- geometric.mean(combined[combined.invariant.set])


##### Red is the NON-ADJUSTED GEOMETRIC MEAN
##### Blue is the Quantiles 2 and 3 ADJUSTED GEOMETRIC MEAN
abline(v=org.gmean, col="red")
abline(v=iqr.gmean, col="blue")
arrows(x0=org.gmean, y0=len+arrow.len, x1=org.gmean, y1=len,length = 0.2, col="red", lwd=4)
arrows(x0=iqr.gmean, y0=len+arrow.len, x1=iqr.gmean, y1=len,length = 0.2, col="blue", lwd=4)
dev.off()

################################### FIGURE S1b #################################
##### Generate the corresponding boxplot of the 150 variates.
##### Blue represents the IQR- geometric mean. Red represents the geometric mean 
##### of the entire set.
################################### FIGURE S1b #################################

bp.name <- paste(figs.dir, "Fig_S1b.png",sep="")
png(bp.name)
boxplot(combined, col="blue")
abline(iqr.gmean,0,col="blue", lwd=2)
abline(org.gmean,0,col="red",lwd=2)
dev.off()