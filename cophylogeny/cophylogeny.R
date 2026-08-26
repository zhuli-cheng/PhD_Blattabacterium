library(ape)
library(phytools)
library(geiger)
library(ggtree)
library(ggplot2)
library(dplyr)
library(dendextend)

source("../TaxonTips/TaxonTips.R")

# input
mito.tree <- read.tree("../MitochondrialTree/input/MT448andRefSeqRY.treefile")
rRNA.tree <- read.tree("./input/RAxML_bestTree.S7A")
protein.tree <- read.tree("./input/244_RY.treefile")
#tRNAsyn.tree <- read.tree("./input/tRNA_synthases_RY.treefile")

# reroot the trees
reroot_tree <- function(tree) {
  two.families <- c("CAM19_1_C019_CR_S19_L001", "YS-3-1_S60_L003")
  common.ancestor <- getMRCA(tree, two.families); common.ancestor
  root(tree, node = common.ancestor, resolve.root = T)
}

mito.tree <- reroot_tree(mito.tree)
rRNA.tree <- reroot_tree(rRNA.tree)
protein.tree <- reroot_tree(protein.tree)
#tRNAsyn.tree <- reroot_tree(tRNAsyn.tree)


# prepare for tanglegram
dendrogram <- function(tree) {
  tree <- chronos(tree, lambda = 1)
  tree <- multi2di(tree)
  tree <- ape::ladderize(tree)
  tree.hclust <- ape::as.hclust.phylo(tree) 
  as.dendrogram(tree.hclust)
}

mito.dend <- dendrogram(mito.tree)
rRNA.dend <-  dendrogram(rRNA.tree)
protein.dend <- dendrogram(protein.tree)
#tRNAsyn.dend <- dendrogram(tRNAsyn.tree)

# tanglegram
#tanglegram(mito.dend, rRNA.dend, 
#           sort = T,
#           highlight_distinct_edges = T, 
#           common_subtrees_color_lines = T, 
#           intersecting = T,
#           main = "Tanglegram of mito.tree and rRNA.tree")


###### congruence between protein tree and rRNA tree
pdf("~/Desktop/Thesis/Typst/figures/chap3results/rRNA_protein_tanglegram.pdf", width = 14, height = 20)
common <- intersect(labels(protein.dend), labels(rRNA.dend))

protein.common <- prune(protein.dend,
                        setdiff(labels(protein.dend), common))

rRNA.common <- prune(rRNA.dend,
                     setdiff(labels(rRNA.dend), common))


dends <- untangle(protein.common, rRNA.common, method = "step2side")

tanglegram(dends[[1]],
          dends[[2]],
           main_left = "Protein-coding genes",
           main_right = "rRNA genes",
           highlight_distinct_edges = TRUE, 
           common_subtrees_color_lines = TRUE, 
           match_order_by_labels = TRUE,
          columns_width = c(8, 3, 8),
           margin_inner = 18,
           lab.cex = 1.2,
           axes = FALSE)
dev.off()

###### congruence between protein tree and mitochondrial tree
pdf("~/Desktop/Thesis/Typst/figures/chap3results/cockroach_Blattabacterium_tanglegram.pdf", width = 14, height = 20)
common <- intersect(labels(protein.dend), labels(mito.dend))

protein.common <- prune(protein.dend,
                        setdiff(labels(protein.dend), common))

mito.common <- prune(mito.dend,
                     setdiff(labels(mito.dend), common))


dends <- untangle(protein.common, mito.common, method = "step2side")

tanglegram(dends[[1]],
           dends[[2]],
           main_left = "Blattabacterium",
           main_right = "Cockroach",
           highlight_distinct_edges = TRUE, 
           common_subtrees_color_lines = TRUE, 
           match_order_by_labels = TRUE,
           columns_width = c(8, 3, 8),
           margin_inner = 18,
           lab.cex = 1.2,
           axes = FALSE)
dev.off()

#plot common tips
common_tips <- intersect(rRNA.tree$tip.label, tRNAsyn.tree$tip.label)
protein.tree.pruned <- drop.tip(protein.tree, setdiff(protein.tree$tip.label, common_tips))
protein.tree.pruned <- ape::ladderize(protein.tree.pruned)
p <- ggtree(protein.tree.pruned) +
  geom_tiplab(size = 1.5, hjust = -0.1) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), legend.position.inside = c(0, 1), 
        legend.justification = c("left", "top"), legend.title = element_blank()); p 
ggsave(filename = "./output/protein_tree_pruned.pdf", plot = p, width = 12, height = 15, dpi = 1200)
system(paste("open", shQuote("./output/protein_tree_pruned.pdf")))

write.tree(protein.tree.pruned, "./output/protein_tree_pruned.treefile")

rRNA.tree.pruned <- drop.tip(rRNA.tree, setdiff(rRNA.tree$tip.label, common_tips))
rRNA.tree.pruned <- ape::ladderize(rRNA.tree.pruned)
p <- ggtree(rRNA.tree.pruned) +
  geom_tiplab(size = 1.5, hjust = -0.1) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), legend.position.inside = c(0, 1), 
        legend.justification = c("left", "top"), legend.title = element_blank()); p 
ggsave(filename = "./output/rRNA_tree_pruned.pdf", plot = p, width = 12, height = 15, dpi = 1200)
system(paste("open", shQuote("./output/rRNA_tree_pruned.pdf")))

tRNAsyn.tree.pruned <- drop.tip(tRNAsyn.tree, setdiff(tRNAsyn.tree$tip.label, common_tips))
tRNAsyn.tree.pruned <- ape::ladderize(tRNAsyn.tree.pruned)
p <- ggtree(tRNAsyn.tree.pruned) +
  geom_tiplab(size = 1.5, hjust = -0.1) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), legend.position.inside = c(0, 1), 
        legend.justification = c("left", "top"), legend.title = element_blank()); p 
ggsave(filename = "./output/tRNAsyn_tree_pruned.pdf", plot = p, width = 12, height = 15, dpi = 1200)
system(paste("open", shQuote("./output/tRNAsyn_tree_pruned.pdf")))

