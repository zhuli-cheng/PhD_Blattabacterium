library(ape)
library(phytools)
library(geiger)
library(ggtree)
library(ggplot2)
library(dplyr)

######### 1. load and trim tree (protein tree, w/o branch length, with highly similar genomes collapsed)
source("Main.R")
source("CollapsingHighlySimilarGenomes.R")

######### 2. define function that modifies edge length (to be replaced with w) based on host family
rename_branches_no_loss <- function(family) {
  tips <- family.presence$tip_label[family.presence$host_family == family]
  ancestral.node <- getMRCA(collapsed.tree, tips)
  descendant_edges <- which(collapsed.tree$edge[, 1] == ancestral.node | collapsed.tree$edge[, 2] %in% getDescendants(collapsed.tree, ancestral.node, "all"))
  collapsed.tree$edge.length[descendant_edges] <<- 1
}

rename_branches_partial_loss <- function(family) {
  tips <- family.presence$tip_label[family.presence$host_family == family]
  ancestral.node <- getMRCA(collapsed.tree, tips)
  descendant_edges <- which(collapsed.tree$edge[, 1] == ancestral.node | collapsed.tree$edge[, 2] %in% getDescendants(collapsed.tree, ancestral.node, "all"))
  collapsed.tree$edge.length[descendant_edges] <<- 2
}

rename_branches_complete_loss <- function(family) {
  tips <- family.presence$tip_label[family.presence$host_family == family]
  ancestral.node <- getMRCA(collapsed.tree, tips)
  descendant_edges <- which(collapsed.tree$edge[, 1] == ancestral.node | collapsed.tree$edge[, 2] %in% getDescendants(collapsed.tree, ancestral.node, "all"))
  collapsed.tree$edge.length[descendant_edges] <<- 3
}

######### 3. for each orthogroup, modify the w 
orthogroups.M2 <- c()

for (i in 1:nrow(orthogroup.sample)) {
  #load gene presence/absence
  orthogroup.name <- paste0("new_orthogroup", i)
  orthogroup.presence <- setNames(ifelse(sample.orthogroup[[orthogroup.name]] == 0, 0, 1), sample.orthogroup$sample)
  
  #find families that have gene loss
  family.presence <- data.frame(tip_label = collapsed.tree$tip.label) %>%
    mutate(host_family = taxonomy$family[match(collapsed.tree$tip.label, taxonomy$sample)]) %>%
    mutate(gene_presence = orthogroup.presence[collapsed.tree$tip.label]) %>%
    group_by(host_family) %>%
    mutate(family_presence = ifelse(sum(gene_presence) == 0, "complete_gene_loss", ifelse(sum(gene_presence) == n(), "no_gene_loss", "partial_gene_loss"))) %>%
    ungroup()
  
  #store orthogroups that have at least one family with partial gene loss and one family with no gene loss
  families.no.loss <- unique(family.presence$host_family[family.presence$family_presence == "no_gene_loss"])
  families.partial.loss <- unique(family.presence$host_family[family.presence$family_presence == "partial_gene_loss"])
  families.complete.loss <- unique(family.presence$host_family[family.presence$family_presence == "complete_gene_loss"])
  if (length(families.partial.loss) > 0 && length(families.no.loss) > 0) {
    orthogroups.M2 <- c(orthogroups.M2, i)
  }
  
  #modify edge length (future w)
  collapsed.tree$edge.length[] <- 0
  lapply(families.no.loss, rename_branches_no_loss)
  lapply(families.partial.loss, rename_branches_partial_loss)
  lapply(families.complete.loss, rename_branches_complete_loss)
  collapsed.tree$node.label <- NULL
  
  #save the tree for bash
  write.tree(collapsed.tree, paste0("../244/output/M2_treefiles/", i, ".tree"))
}

cat(orthogroups.M2, sep = "\n")

######### analyse CodeML results
dN.tree <- read.tree("../244/input/M2_1_dN.treefile")
plot(dN.tree)
dS.tree <- read.tree("../244/input/M2_1_dS.treefile")
plot(dS.tree)
w.tree <- read.tree("../244/input/M2_1_w.treefile")
plot(w.tree)

summary.M2 <- read.csv("../244/input/M2_by_gene_loss_summary.csv", header = F) %>%
  setNames(c("new_orthogroup", "w0", "w1")) %>%
  left_join(orthogroup.sample)

ggplot(summary.M2, aes(x = w0, y = w1)) +
  geom_point(size = 1) +
  facet_wrap(~ curated_COG_category) +
  labs(
    title = "Comparison of w0 and w1 by curated_COG_category",
    x = "w0",
    y = "w1"
  ) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") # Add y = x line


