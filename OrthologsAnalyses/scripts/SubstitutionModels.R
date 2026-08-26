library(ape)
library(phytools)
library(geiger)
library(ggtree)
library(ggplot2)
library(dplyr)

source("../../TaxonTips/TaxonTips.R")

####################### read in input data #################################
sample.size <- 244

#read in phylogeny 
GTR.tree <- read.tree(paste0("../", sample.size, "plus9/input/244_GTR.treefile"))
aa.tree <- read.tree(paste0("../", sample.size, "plus9/input/244plus9_aa.treefile"))
nonRev.tree <- read.tree(paste0("../", sample.size, "plus9/input/244plus9_nonReversible.treefile"))

#reroot the tree, using outgroups
reroot <- function(tree){
  two.families <- c("GCF_000210835.1", "GCF_000236705.1")
  common.ancestor <- getMRCA(tree, two.families); common.ancestor
  rerooted.tree <<- root(tree, node = common.ancestor, resolve.root = T)
}

reroot(GTR.tree)
reroot(aa.tree)
reroot(nonRev.tree)

#collapse highly similar genomes; run when needed
collapsed.tree <- rerooted.tree
source("CollapsingHighlySimilarGenomes.R")

############################ visualize the tree #############################
#using host family
modify_tip(collapsed.tree, "family")
family.tip.tree <- ggtree(collapsed.tree) %<+% taxonomy.tip +
  geom_tiplab(size = 2, hjust = -0.1, aes(label = taxon)); family.tip.tree
ggsave(filename = "../244plus9/output/family_tip_treeaa.pdf", plot = family.tip.tree, width = 10, height = 15, dpi = 1200)

plot.phylo(collapsed.tree, show.node.label = TRUE, edge.color = "blue")

#using host taxa
modify_tip(collapsed.tree, "mixed")
mixed.tip.tree <- ggtree(collapsed.tree) %<+% taxonomy_tip +
  geom_tiplab(size = 2, hjust = -0.1, aes(label = taxon)); mixed.tip.tree
ggsave(filename = "../244plus9/output/mixed_tip_tree.pdf", plot = mixed.tip.tree, width = 10, height = 15, dpi = 1200)

#using sample id
ID.tip.tree <- ggtree(collapsed.tree) +
  geom_tiplab(size = 2, hjust = -0.1); ID.tip.tree
ggsave(filename = "../244plus9/output/ID_tip_tree.pdf", plot = ID.tip.tree, width = 10, height = 15, dpi = 1200)


