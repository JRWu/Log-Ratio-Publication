# Figure 5 with MGNify MGYS00001863 dataset

# https://www.ebi.ac.uk/metagenomics/projects/ERP023694
# https://www.ebi.ac.uk/metagenomics/studies/MGYS00001863 for actual dataset
# https://www.biorxiv.org/content/biorxiv/early/2018/05/11/248302.full.pdf

eff.plot <- function(x, main=""){
	plot(x$diff.win, x$diff.btw, pch=19, col=rgb(0.8,0.7,0.5,0.3), cex=0.5,
	xlim=c(0,8), ylim=c(-7,10), xlab="Dispersion", ylab="Difference", main=main)
	points(x[sparse.set,"diff.win"], x[sparse.set,"diff.btw"], pch=19, col=rgb(0.2,0.1,0.05,0.3), cex=0.5)
	points(x[ribo,"diff.win"], x[ribo,"diff.btw"], pch=19, col=rgb(0,0,1,1), cex=1)
	points(x[glycol,"diff.win"], x[glycol,"diff.btw"], pch=19, col=rgb(1,0,0.8,1), cex=1)
	abline(0,1, lty=2, lwd=3, col=rgb(0,0,0,0.4))
	abline(0,-1,lty=2, lwd=3, col=rgb(0,0,0,0.4))
	abline(h=0,lty=3, lwd=3, col=rgb(0,0,0,0.4))
	text(6,9.5, labels=paste("transl=", round(median(x[ribo,"diff.btw"], na.rm=T),3), sep=""))
	text(6,8.5, labels=paste("tRNA=", round(median(x[glycol,"diff.btw"],na.rm=T),3), sep=""))
}

library(CoDaSeq)
library(ALDEx2)
d <- read.table("data/ERP023694_GO_abundances_v3.0.tsv", header=T, row.names=1, sep="\t", quote="")
m <- read.table("data/pub.metadata.txt", header=T, row.names=1, sep="\t", comment.char="", quote="")

d.annot <- data.frame(d[,1])
d.category <- data.frame(d[,2])
rownames(d.annot) <- rownames(d)
rownames(d.category) <- rownames(d)

d[,1:2] <- NULL

# keep only biological process GO annotations
d.proc <- d[grep("process", d.category[,1]),]

# choose H and BV samples
d.var <- codaSeq.filter(d.proc, var.filt=TRUE, samples.by.row=FALSE)
var.clr <- apply(d.var +0.5, 2, function(x) log(x) - mean(log(x)))
var.pcx <- prcomp(t(var.clr))
H <- which(var.pcx$x[,1] < -10)
B <- which(var.pcx$x[,1] > 10)

d.ord <- data.frame(d.proc[,B], d.proc[,H])
sparse.set <- rownames(d.ord)[apply(d.ord, 1, min) == 0]

conds <- c(rep("B", length(B)), rep("H", length(H)))

user <- which(rownames(d.ord) %in% ribo)

aldex_plot <- function(denom, main=""){
    x <- aldex.clr(d.ord, conds, denom=denom)
    x.e <- aldex.effect(x)
    transl <- rownames(d.annot)[grep("translat", d.annot[,1])]
    ribo <- transl[x.e[transl, "rab.all"] >  quantile(x.e$rab.all)[4]]
    tRNA <- rownames(d.annot)[grep("tRNA", d.annot[,1])]
    glycol <- tRNA[x.e[tRNA, "rab.all"] >  quantile(x.e$rab.all)[4]]

    eff.plot(x.e, main=main)
}
pdf("figures/mgnify.pdf", height=10, width=9)
par(mfrow=c(2,2))
aldex_plot("all", "all")
aldex_plot("lvha", "lvha")
aldex_plot("iqlr", "iqlr")
aldex_plot(user, "transl")
dev.off()
