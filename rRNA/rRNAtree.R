library(ape)
library(phytools)
library(geiger)
library(ggtree)
library(ggplot2)
library(dplyr)

#input
rRNA.tree <- read.tree("../rRNA/input/RAxML_bestTree.S7A")
circularisation <- read.csv("../rRNA/input/circularisation.csv")

reroot_tree <- function(tree) {
  two.families <- c("CAM19_1_C019_CR_S19_L001", "YS-3-1_S60_L003")
  common.ancestor <- getMRCA(tree, two.families); common.ancestor
  root(tree, node = common.ancestor, resolve.root = T)
}

rRNA.tree <- reroot_tree(rRNA.tree)


p <- ggtree(rRNA.tree) +
  geom_tiplab(size = 1.5, hjust = -0.1); p 

#modify the tip labels with family information
source("../TaxonTips/TaxonTips.R")
modify_tip(rRNA.tree, "mixed")

p <- ggtree(rRNA.tree, ladderize = FALSE) %<+% taxonomy.tip +
  geom_tiplab(size = 3, hjust = -0.1, aes(label = taxon)) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +
  geom_treescale(x = 0, y = 0) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.position.inside = c(0, 1), 
        legend.justification = c("left", "top"), legend.title = element_blank()); p #no node support
ggsave(filename = "../rRNA/output/rRNA_tree.pdf", plot = p, width = 10, height = 10, dpi = 1200)
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/Blattabacterium_rRNAs_mixed_tip_tree.pdf", plot = p, width = 10, height = 14, dpi = 1200)

system(paste("open", shQuote("../rRNA/output/rRNA_tree.pdf")))


#add circularisation information on the tips
taxonomy.tip.circularisation <- taxonomy.tip %>% 
  left_join(circularisation, by = c("tip_label" = "Sample")) %>%
  select(tip_label, taxon, circular) %>%
  mutate(condition = ifelse(grepl("#N/A", circular), "uncircularised", "circular"))

p <- ggtree(rRNA.tree) %<+% taxonomy.tip.circularisation +
  geom_tiplab(size = 1.5, hjust = -0.1, aes(label = taxon, color = condition)) +
  scale_color_manual(values = c("uncircularised" = "red", "circular" = "black")) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.position.inside = c(0, 1), 
        legend.justification = c("left", "top"), legend.title = element_blank()); p 
#no interesting signal
ggsave(filename = "../rRNA/output/tree_with_circularisation.pdf", plot = p, width = 10, height = 16, dpi = 1200)

system(paste0("open ", "../rRNA/output/tree_with_circularisation.pdf"))
