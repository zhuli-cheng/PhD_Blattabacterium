source("Main.R")
#make sure the "AncestralStateReconstruction.R" has been commented out, to maintain consistency within clades

################# show each sample's number of orthogroups on the tree ###############
#specify orthogroups to look at
no.of.orthogroups <- setNames(as.numeric(sample.orthogroup$no_of_orthologs), sample.orthogroup$sample)
all(names(no.of.orthogroups) %in% tree$tip.label)
all(tree$tip.label %in% names(no.of.orthogroups))


#map statistics on the tree tips
modify_tip(collapsed.tree, "mixed")
taxonomy.vector <- setNames(taxonomy$mixed, taxonomy$sample)
collapsed.tree$tip.label <- paste0(taxonomy.vector[collapsed.tree$tip.label], ":", no.of.orthogroups[collapsed.tree$tip.label])
tip.data <- data.frame(tip_label = collapsed.tree$tip.label)

no.of.orthogroups.GC <- ggtree(collapsed.tree, ladderize = FALSE) %<+% tip.data +
  geom_tiplab(size = 3, hjust = -0.1) +
  geom_treescale(x = 0, y = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.position = c(0, 1), legend.justification = c("left", "top"), legend.title = element_blank()); no.of.orthogroups.GC

ggsave("~/Desktop/Thesis/Typst/figures/chap3results/genome_gene_loss.png", plot = no.of.orthogroups.GC, width = 10, height = 10, dpi = 1200)

#deletions on YS-3
YS <- orthogroup.sample %>%
  select("gene", "YS-3-1_S60_L003", "YS-5_S62_L003") %>%
  filter(if_any(everything(), is.na))

         