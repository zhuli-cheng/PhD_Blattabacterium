source("Main.R")

#full tree
ID.tip.tree <- ggtree(mito.tree) %<+% absence +
  geom_tiplab(size = 2.5, hjust = -0.1, aes(color = absence_confirmed)) +  
  geom_treescale(x = 0, y = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +
  geom_text2(aes(subset = !isTip, label = label),
             hjust = -0.3,
             vjust = -0.3,
             size = 2); ID.tip.tree
ggsave(filename = "../output/Blattabacterium_ID_tip_tree.pdf", plot = ID.tip.tree, width = 20, height = 50, dpi = 1200, limitsize = FALSE)


#collapse highly similar genomes
source("./CollapsingHighlySimilarGenomes.R")
modify_tip(mito.tree, "mixed")

absence <- absence %>% 
  full_join(taxonomy.tip, join_by(sample == tip_label)) %>%
  mutate(taxon = ifelse(is.na(taxon) | taxon == "", sample, taxon))

ID.tip.tree <- ggtree(mito.tree) %<+% absence +
  geom_tiplab(size = 2, hjust = -0.1, aes(label = taxon, color = absence_confirmed)); ID.tip.tree
ggsave(filename = "../output/Wolbachia_ID_tip_tree.pdf", plot = ID.tip.tree, width = 10, height = 16, dpi = 1200)


########## map Wolbachia and Rickettsia data
symbionts <- GC.absence %>%
  select(blattabacterium, wolbachia, rickettsia)
symbionts <- symbionts[mito.tree$tip.label, , drop = FALSE]

phylo.heatmap(mito.tree, symbionts, fsize = 0.01, standardize = TRUE, legend = FALSE, labels = FALSE, split=c(0.8, 0.2)) +
  #add.scale.bar(cex = 0.1, pos = "bottomleft")
  str(symbionts)

