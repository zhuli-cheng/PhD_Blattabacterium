stree <- drop.tip(birthdeath.tree(b = 1, d = 0, taxa.stop = 51), "51")
plot.phylo(stree)
plotTree(stree)
stree <- rescaleTree(stree, 30)
x <- fastBM(stree, sig2 = 0.1, bounds = c(0, 1), internal = T); x
phenogram(stree, x[1:51], ftype = "reg")

Q <- matrix(c(-1, 1, 10, -10), 2, 2)
mtree <- sim.history(stree, Q)
cols <- c("blue", "pink")
names(cols) <- c (1, 2)
plotSimmap(mtree, cols)
nchanges <- sum(sapply(mtree$maps, length)) - nrow(stree$edge); nchanges

sig2 <- c(10, 10)
names(sig2) <- c(1, 2)
x <- sim.rates(mtree, sig2)
layout(c(1, 2))
plotSimmap(mtree)
phenogram(stree, x)

tips <- mtree$states; tips