

#map tip colors
ordered.trait.data <- orthogroup426[rerooted.tree$tip.label]
colors <- ifelse(ordered.trait.data == 0, "black", "green")
tip.colors <- colors[as.character(ordered.trait.data)]
plot(rerooted.tree, show.tip.label = T, tip.color = colors)
#write.tree(rerooted.tree, file = "tree.nwk")


#reconstruct ancestral state

orthogroup426binary <- ifelse(orthogroup426 == 0, 0, 1)
simmap <- make.simmap(rerooted.tree, orthogroup426binary, model = Q, nsim = 1)
fit <- fitDiscrete(rerooted.tree, orthogroup426binary, model = Q)
summary <- describe.simmap(simmap)
plotSimmap(simmap, fsize = 0.7)
#tiplabels(pch = 21, bg = tip.colors, cex = 0.7)

#write.tree(simmap$tree, file = "simulated gene loss.nwk")

p <- ggtree(collapsed.tree) +
  geom_text2(aes(label = node), hjust = -0.3)
viewClade(p, 329)

library(ape)
library(phytools)
library(geiger)

tree <- read.tree("245_GTR.treefile")
gene <- read.csv("sample.orthogroup.csv", header = T)

gene24 <- setNames(as.character(gene$V24), gene$X)
gene500 <- setNames(as.character(gene$V500), gene$X)


#reroot the tree
tip_of_interest <- c("CAM19_2_C166_CR_S65_L004_Corydiidae", "CAM19_1_C019_CR_S19_L001_Anaplectidae")
common_ancestor <- getMRCA(tree, tip_of_interest); common_ancestor
rerooted_Tree <- root(tree, node = common_ancestor, resolve.root = TRUE)
plot(rerooted_Tree)
plot.phylo(tree)
write.tree(rerooted_Tree, file="tree.nwk")

#remove some samples
rerooted_Tree <- drop.tip(rerooted_Tree, "CAM19_1_C074_CR_S73_L001_undetermined")
all(names(gene24) %in% rerooted_Tree$tip.label)
all(rerooted_Tree$tip.label %in% names(gene24))

#specify a matrix where 0 cannot transition to 1
Q <- matrix(c(-1, 0, 1, -1), nrow = 2)
rownames(Q) <- colnames(Q) <- c("0", "1")

#map tip colors
ordered_trait_data <- gene24[tree$tip.label]
colors <- c("0" = "black", "1" = "green")
tip_colors <- colors[as.character(ordered_trait_data)]
plot(rerooted_Tree, tip.color = tip_colors, main = "Tip States Based on Trait Data")



#reconstruct ancestral state
result <- ace(gene24, rerooted_Tree, model = Q)
simmap <- make.simmap(rerooted_Tree, gene24, model = Q, nsim = 1)
fit <- fitDiscrete(rerooted_Tree, gene24, model = Q)
summary <- describe.simmap(simmap)
plotSimmap(simmap, fsize = 0.7, colors = colors)
tiplabels(pch = 21, bg = tip_colors, cex = 0.7)

#collapse short branches
threshold <- 0.01

collapse_short_branches <- function(tree, threshold) {
  branch_lengths <- tree$edge.length
  short_branches <- which(branch_lengths < threshold)
  if (length(short_branches) > 0) {
    tree <- drop.tip(tree, which(tree$edge.length < threshold))
  }
  return(tree)
}

collapsed_tree <- collapse_short_branches(rerooted_Tree, threshold)

plot(collapsed_tree, main = "Tree with Short Branches Collapsed")


