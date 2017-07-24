# see 0_git/twntyfr/chunk/ALDEx_corrections.rablibrary(zCompositions)
library(ALDEx2)

e.min <- read.table("data/twntyfr.txt", header=T,    row.names=1, check.names=F, sep="\t", comment.char="", quote="")

ribo <- c(grep("LSU", rownames(e.min)), grep("SSU", rownames(e.min)))
glycol <- c(2418,1392,1305,1306,2421,1049)
sparse.set <- names(which(apply(e.min, 1, min) == 0))

e.min.n0.CZM <- cmultRepl(t(e.min), label=0, method="CZM")
e.min.clr <- t( apply(e.min.n0.CZM, 1, function(x){log(x) - mean(log(x))}) )

e.min.pcx <- prcomp(e.min.clr)

e.min.clr.u <- t( apply(e.min.n0.CZM, 1, function(x){log(x) - mean(log(x[ribo]))}) )

e.min.pcx.u <- prcomp(e.min.clr.u)

conds <-c("H","H","H","H","B","H","B","B","H","B","H","B","B","B","B","B","B","B","H","B","H","H")


eff.plot <- function(x, main=""){
	plot(x$diff.win, x$diff.btw, pch=19, col=rgb(0.8,0.7,0.5,0.3), cex=0.5, xlab="Dispersion", ylab="Difference", main=main)
	points(x[sparse.set,"diff.win"], x[sparse.set,"diff.btw"], pch=19, col=rgb(0.2,0.1,0.05,0.3), cex=0.5)
	points(x$diff.win[ribo], x$diff.btw[ribo], pch=19, col=rgb(0,0,1,1), cex=0.5)
	points(x$diff.win[glycol], x$diff.btw[glycol], pch=19, col=rgb(1,0,0.8,1), cex=0.5)
	abline(0,1, lty=2, lwd=3, col=rgb(0,0,0,0.4))
	abline(0,-1,lty=2, lwd=3, col=rgb(0,0,0,0.4))
	abline(h=0,lty=3, lwd=3, col=rgb(0,0,0,0.4))
}

ma.plot <- function(x, x1, y1,y2, main=""){
	plot(x$rab.all, x$diff.btw, pch=19, col=rgb(0.8,0.7,0.5,0.3), cex=0.5, xlab="Log-ratio abundance", ylab="Difference", main=main)
	points(x[sparse.set,"rab.all"], x[sparse.set,"diff.btw"], pch=19, col=rgb(0.2,0.1,0.05,0.3), cex=0.5)
	points(x$rab.all[ribo], x$diff.btw[ribo], pch=19, col=rgb(0,0,1,1), cex=0.5)
	points(x$rab.all[glycol], x$diff.btw[glycol], pch=19, col=rgb(1,0,0.8,1), cex=0.5)
	abline(h=0,lty=3, lwd=3, col=rgb(0,0,0,0.4))
	text(x1,y1, labels=paste("med ribo=", round(median(x$diff.btw[ribo]),3), sep=""))
	text(x1,y2, labels=paste("med glyc=", round(median(x$diff.btw[glycol]),3), sep=""))
}


pdf("figures/twtyfr.pdf", height=5,width=9)
par(mfrow=c(1,2))
#x <- aldex.clr(e.min, conds)
#x.e <- aldex.effect(x, conds, verbose=FALSE)
eff.plot(x.e, main="Effect plot")
legend(6,-6, legend=c("non-sparse", "sparse", "ribosome", "glycolysis"),   col=c(rgb(0.8,0.7,0.5,0.3), rgb(0.2,0.1,0.05,0.3), rgb(0,0,1,1), rgb(1,0,0.8,1)), pch=19, cex=0.5, bg="white")

ma.plot(x.e, 10,9,7.5, main="Bland Altman plot")

dev.off()

pdf("figures/twtyfr_corr.pdf", height=9,width=12)
par(mfrow=c(2,3))

#x <- aldex.clr(e.min, conds, denom="iqlr")
#x.i <- aldex.effect(x, conds, verbose=FALSE)
eff.plot(x.i, main="Effect plot: IQLR")
	legend(6,-6, legend=c("non-sparse", "sparse", "ribosome", "glycolysis"),   col=c(rgb(0.8,0.7,0.5,0.3), rgb(0.2,0.1,0.05,0.3), rgb(0,0,1,1), rgb(1,0,0.8,1)), pch=19, cex=1, bg="white")

#x <- aldex.clr(e.min, conds, denom="zero")
#x.z <- aldex.effect(x, conds, verbose=FALSE)
eff.plot(x.z, main="Effect plot: Zero")

#x <- aldex.clr(e.min, conds, denom=ribo)
#x.u <- aldex.effect(x, conds, verbose=FALSE)
eff.plot(x.u, main="Effect plot: Ribosomal functions")

#x <- aldex.clr(e.min, conds, denom="iqlr")
#x.i <- aldex.effect(x, conds, verbose=FALSE)
ma.plot(x.i, 10,9,7.5,main="BA plot: IQLR")

#x <- aldex.clr(e.min, conds, denom="zero")
#x.z <- aldex.effect(x, conds, verbose=FALSE)
ma.plot(x.z, 10,9,7.5,main="BA plot: Zero")

#x <- aldex.clr(e.min, conds, denom=ribo)
#x.u <- aldex.effect(x, conds, verbose=FALSE)
ma.plot(x.u, -0,9,7.5,main="BA plot: Ribosomal functions")

dev.off()


	text(x1,y1, labels=paste("med ribo=", round(median(x$diff.btw[ribo]),3), sep=""))
	text(x1,y2, labels=paste("med glyc=", round(median(x$diff.btw[glycol]),3), sep=""))


pdf("figures/twtyfr_poster.pdf", height=3.7,width=12)
par(mfrow=c(1,4))

eff.plot(x.e, main="Denom=all")
legend(5.6,-4.1, legend=c("abundant", "rare", "ribosome", "glycolysis"),   col=c(rgb(0.8,0.7,0.5,0.3), rgb(0.2,0.1,0.05,0.3), rgb(0,0,1,1), rgb(1,0,0.8,1)), pch=19, cex=1, bg="white")
	text(5.5,10, labels=paste("med ribo=", round(median(x.e$diff.btw[ribo]),3), sep=""))
	text(5.5,9, labels=paste("med glyc=", round(median(x.e$diff.btw[glycol]),3), sep=""))

#x <- aldex.clr(e.min, conds, denom="iqlr")
#x.i <- aldex.effect(x, conds, verbose=FALSE)
eff.plot(x.i, main="Denom=IQLR")
#	legend(6,-6, legend=c("non-sparse", "sparse", "ribosome", "glycolysis"),   col=c(rgb(0.8,0.7,0.5,0.3), rgb(0.2,0.1,0.05,0.3), rgb(0,0,1,1), rgb(1,0,0.8,1)), pch=19, cex=1, bg="white")
	text(5.5,9.5, labels=paste("med ribo=", round(median(x.i$diff.btw[ribo]),3), sep=""))
	text(5.5,8.5, labels=paste("med glyc=", round(median(x.i$diff.btw[glycol]),3), sep=""))

#x <- aldex.clr(e.min, conds, denom="zero")
#x.z <- aldex.effect(x, conds, verbose=FALSE)
eff.plot(x.z, main="Denom=non-zero")
	text(5.5,9.5, labels=paste("med ribo=", round(median(x.z$diff.btw[ribo]),3), sep=""))
	text(5.5,8.5, labels=paste("med glyc=", round(median(x.z$diff.btw[glycol]),3), sep=""))

#x <- aldex.clr(e.min, conds, denom=ribo)
#x.u <- aldex.effect(x, conds, verbose=FALSE)
eff.plot(x.u, main="Denom=Ribosomal functions")
	text(5.5,9.5, labels=paste("med ribo=", round(median(x.u$diff.btw[ribo]),3), sep=""))
	text(5.5,8.5, labels=paste("med glyc=", round(median(x.u$diff.btw[glycol]),3), sep=""))

dev.off()
