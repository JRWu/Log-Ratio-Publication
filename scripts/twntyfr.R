# see 0_git/twntyfr/chunk/ALDEx_corrections.rablibrary(zCompositions)
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

x <- aldex.clr(e.min, conds)
x.e <- aldex.effect(x)

x.l <- aldex.clr(e.min, conds, denom="lvha")
x.e.l <- aldex.effect(x.l)

x.i <- aldex.clr(e.min, conds, denom="iqlr")
x.e.i <- aldex.effect(x.i)

x.r <- aldex.clr(e.min, conds, denom=ribo)
x.e.r <- aldex.effect(x.r)

eff.plot <- function(x, main=""){
	plot(x$diff.win, x$diff.btw, pch=19, col=rgb(0.8,0.7,0.5,0.3), cex=0.5, xlab="Dispersion", ylab="Difference", main=main)
	points(x[sparse.set,"diff.win"], x[sparse.set,"diff.btw"], pch=19, col=rgb(0.2,0.1,0.05,0.3), cex=0.5)
	points(x$diff.win[ribo], x$diff.btw[ribo], pch=19, col=rgb(0,0,1,1), cex=0.5)
	points(x$diff.win[glycol], x$diff.btw[glycol], pch=19, col=rgb(1,0,0.8,1), cex=0.5)
	abline(0,1, lty=2, lwd=3, col=rgb(0,0,0,0.4))
	abline(0,-1,lty=2, lwd=3, col=rgb(0,0,0,0.4))
	abline(h=0,lty=3, lwd=3, col=rgb(0,0,0,0.4))
	text(2,9, labels=paste("ribo=", round(median(x$diff.btw[ribo]),3), sep=""))
	text(2,7, labels=paste("glyc=", round(median(x$diff.btw[glycol]),3), sep=""))
}

ma.plot <- function(x, x1, denom, y1,y2, main=""){
	plot(x$rab.all, x$diff.btw, pch=19, col=rgb(0.8,0.7,0.5,0.3), cex=0.5, xlab="Log-ratio abundance", ylab="Difference", main=main)
	points(x[sparse.set,"rab.all"], x[sparse.set,"diff.btw"], pch=19, col=rgb(0.2,0.1,0.05,0.3), cex=0.5)
	points(x$rab.all[ribo], x$diff.btw[ribo], pch=19, col=rgb(0,0,1,1), cex=0.5)
	points(x$rab.all[glycol], x$diff.btw[glycol], pch=19, col=rgb(1,0,0.8,1), cex=0.5)
	abline(h=0,lty=3, lwd=3, col=rgb(0,0,0,0.4))
	text(x1,y1, labels=paste("med ribo=", round(median(x$diff.btw[ribo]),3), sep=""))
	text(x1,y2, labels=paste("med glyc=", round(median(x$diff.btw[glycol]),3), sep=""))
}

pdf("figures/MAx4.pdf", height=10, width=9)
par(mfrow=c(2,2))
eff.plot(x.e, main="denom = all")
eff.plot(x.e.i, main="denom = IQLR")
eff.plot(x.e.l, main="denom = LVHA")
eff.plot(x.e.r,  main="denom = ribosome")
dev.off()

