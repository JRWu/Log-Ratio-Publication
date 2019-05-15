library(ALDEx2)

e.min <- read.table("data/twntyfr.txt", header=T,    row.names=1, check.names=F, sep="\t", comment.char="", quote="")

ribo <- c(grep("LSU", rownames(e.min)), grep("SSU", rownames(e.min)))
glycol <- c(2418,1392,1305,1306,2421,1049)
sparse.set <- names(which(apply(e.min, 1, min) == 0))

#e.min.n0.CZM <- cmultRepl(t(e.min), label=0, method="CZM")
#e.min.clr <- t( apply(e.min.n0.CZM, 1, function(x){log(x) - mean(log(x))}) )
#
#e.min.pcx <- prcomp(e.min.clr)
#
#e.min.clr.u <- t( apply(e.min.n0.CZM, 1, function(x){log(x) - mean(log(x[ribo]))}) )
#
#e.min.pcx.u <- prcomp(e.min.clr.u)

conds <-c("H","H","H","H","B","H","B","B","H","B","H","B","B","B","B","B","B","B","H","B","H","H")
x <- aldex.clr(e.min, conds=conds, denom="all")
x.e <- aldex.effect(x, conds)

pdf("ngs_BA_hist.pdf", width=8, height=5)
par(mfrow=c(1,2))
plot(x.e$rab.all, x.e$diff.btw, pch=19, col=rgb(0.8,0.7,0.5,0.3), cex=0.5, xlab="Log-ratio abundance", ylab="Difference", main="Bland-Altman")
points(x.e[sparse.set,"rab.all"], x.e[sparse.set,"diff.btw"], pch=19, col=rgb(0.2,0.1,0.05,0.3), cex=0.5)
legend(6, 11, legend=c("abundant", "rare"), pch=19, col=c(rgb(0.8,0.7,0.5,0.3),rgb(0.2,0.1,0.05,0.3)))

hist(x.e$diff.btw, breaks=49, main="Log2 difference", col=rgb(0.8,0.7,0.5,0.3), xlab="Bin value")

dev.off()
ribo <- c(grep("LSU", rownames(e.min)), grep("SSU", rownames(e.min)))
glycol <- c(2418,1392,1305,1306,2421,1049)
sparse.set <- names(which(apply(e.min, 1, min) == 0))

ma.plot <- function(x, x1, y1,y2, main=""){
	plot(x$rab.all, x$diff.btw, pch=19, col=rgb(0.8,0.7,0.5,0.3), cex=0.5, xlab="Log-ratio abundance", ylab="Difference", main=main)
	points(x[sparse.set,"rab.all"], x[sparse.set,"diff.btw"], pch=19, col=rgb(0.2,0.1,0.05,0.3), cex=0.5)
	points(x$rab.all[ribo], x$diff.btw[ribo], pch=19, col=rgb(0,0,1,1), cex=0.5)
	points(x$rab.all[glycol], x$diff.btw[glycol], pch=19, col=rgb(1,0,0.8,1), cex=0.5)
	abline(h=0,lty=3, lwd=3, col=rgb(0,0,0,0.4))
	text(x1,y1, labels=paste("med ribo=", round(median(x$diff.btw[ribo]),3), sep=""))
	text(x1,y2, labels=paste("med glyc=", round(median(x$diff.btw[glycol]),3), sep=""))
}
pdf("ngs_BA_uncentered.pdf", width=7, height=7)
ma.plot(x.e, main="Non-centred")
dev.off()
####### edgeR
library(edgeR)
y <- DGEList(counts=e.min,group=factor(conds))

y <- calcNormFactors(y)

y <- estimateDisp(y)

et <- exactTest(y)
# hack to get the data for plotting later
tt <- topTags(et, sort.by="none", n=6236)

summary(et <- decideTestsDGE(et))
detags <- rownames(y)[as.logical(et)]


# you can plot the MA plot

# compare to the MA plot from ALDEx, they are rather similar for this dataset
pdf("ngs_edgeR.pdf", width=6, height=6)
plotSmear(y, de.tags=detags, main="edgeR",  lowess=T)
abline(h=c(-1, 1), col="blue", lty=2)
dev.off()

