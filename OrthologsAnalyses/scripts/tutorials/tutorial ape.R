library(ggtree)
library(ape)
library(phytools)



p <- ggtree(tree)
p + xlim_tree(0, 30)  # Adjust xlim_tree for zooming in horizontally

subtree <- extract.clade(tree, node = common_ancestor)  # Replace 20 with the node of interest
plot(tree, type="fan")


subtree <- drop.tip(tree, setdiff(tree$tip.label, tip_of_interest))
plot(subtree)

plot(tree)
nodelabels()

tip_of_interest <- c("CAM19_2_C100_CR_S4_L004_Anaplectidae", "CAM19_2_C166_CR_S65_L004_Corydiidae", "AU013_S342_Pseudophyllodromiidae", "CAM19_2_C102_CR_S6_L004_Blaberidae")
subtree <- drop.tip(tree, tip = tip_of_interest, trim.internal = TRUE, subtree = FALSE, collapse.singles = TRUE)
plot(subtree)
nodelabels()



clade <- extract.clade(tree, tip = tip_of_interest)


rerootTree <- root(myTree, node = 7)
plot(rerootTree, node.depth = 2)

nodelabels()

edgecol <- c()
edgecol[1:20] = "red"

tipcol <- c()

rerootTree$tip.label <- c()

plot.phylo(rerootTree, "fan", 
           use.edge.length = TRUE,
           main = "title", 
           edge.color = edgecol,
           tip.color = tipcol,
           show.tip.label = TRUE,
           font = 3,
           cex = 1)
