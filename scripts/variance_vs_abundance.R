# plots of relative abundance vs. dispersion in different datasets

library(ALDEx2)
# selex
data(selex)
conds <- c(rep("A",7), rep("B",7))
x <- aldex.clr(selex,conds)
x.lvha <- aldex.clr(selex, conds, denom="lvha")
x.iqlr <- aldex.clr(selex, conds, denom="iqlr")
x.zero <- aldex.clr(selex, conds, denom="zero")

x.e <- aldex.efffect(x)

 t4 <-  read.table("/Users/ggloor/Documents/0_git/Log-Ratio-Publication/data/twntyfr.txt", header=T, sep="\t", quote="", comment.char="", row.names=1)

conds.t <-c("H","H","H","H","B","H","B","B","H","B","H","B","B","B","B","B","B","B","H","B","H","H")

t <- aldex.clr(t4, conds.t)
t.e <- aldex.effect(t)
t.lvha <- aldex.clr(t4, conds.t, denom="lvha")
t.zero <- aldex.clr(t4, conds.t, denom="zero")
t.iqlr <- aldex.clr(t4, conds.t, denom="iqlr")

# tiyani
load("/Users/ggloor/Documents/0_git/effect/data/ref.set.ty.Rdata")

# barton yeast
load("/Users/ggloor/Documents/0_git/effect/data/ref.set.yeast.Rdata")

pdf("figures/Disp-v-Abund.pdf")
par(mfrow=c(1,1))

plot(ref.set.ty[[1]]$diff.win, ref.set.ty[[1]]$rab.all,col=rgb(1,0,0,0.2), pch=19,
   xlim=c(0,6), ylim=c(-10,10), xlab="Dispersion", ylab="log2 Relative Abundance")

points(ref.set.yeast[[1]]$diff.win, ref.set.yeast[[1]]$rab.all, col=rgb(0,0,1,0.2), pch=19)
legend(0,-7, legend=c("transcriptome", "microbiome"), col=c("blue","red"), pch=19)
dev.off()

plot(x.e$diff.win, x.e$rab.all, col=rgb(0,0,0,0.2), pch=19)
points(x.e$diff.win[x.iqlr@denom], x.e$rab.all[x.iqlr@denom], col=rgb(1,0,0,0.5), pch=19)


plot(t.e$diff.win, t.e$rab.all, col=rgb(0,0,0,0.2),pch=19)
points(t.e$diff.win[t.iqlr@denom], t.e$rab.all[t.iqlr@denom], col=rgb(1,0,0,0.5), pch=19)
points(t.e$diff.win[t.lvha@denom], t.e$rab.all[t.lvha@denom], col=rgb(1,0,1,0.5), pch=19)
points(t.e$diff.win[t.zero@denom[[1]]], t.e$rab.all[t.zero@denom[[1]]], col=rgb(0,1,1,0.5), pch=19)




