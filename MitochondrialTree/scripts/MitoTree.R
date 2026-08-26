source("./Main.R")

#using host family
modify_tip(mito.tree, "family")
family.tip.tree <- ggtree(mito.tree) %<+% taxonomy.tip +
  geom_tiplab(size = 1.5, hjust = -0.1, aes(label = taxon)) +
  geom_treescale(x = 0, y = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +
  geom_text2(aes(subset = !isTip, label = label),
             hjust = -0.3,
             vjust = -0.3,
             size = 2); family.tip.tree
ggsave(filename = "../output/family_tip_tree.pdf", plot = family.tip.tree, width = 10, height = 15, dpi = 1200)

#using host taxa
modify_tip(mito.tree, "mixed")
mixed.tip.tree <- ggtree(mito.tree) %<+% taxonomy.tip +
  geom_tiplab(size = 2, hjust = -0.1, aes(label = taxon)) +
  geom_treescale(x = 0, y = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +
  geom_text2(aes(subset = !isTip, label = label),
             hjust = -0.3,
             vjust = -0.3,
             size = 2); mixed.tip.tree
ggsave(filename = "../output/mixed_tip_tree.pdf", plot = mixed.tip.tree, width = 12, height = 16, dpi = 1200)

#collapse tree for thesis
modify_tip(mito.tree, "mixed")
d <- cophenetic.phylo(mito.tree)
hc <- hclust(as.dist(d), method = "complete")
groups <- cutree(hc, h = 0.05)
rep_tips <- tapply(names(groups), groups, `[`, 1)
collapsed.tree <- drop.tip(mito.tree, setdiff(mito.tree$tip.label, rep_tips))

mix.tip.tree.collapsed <- ggtree(collapsed.tree, ladderize = FALSE) %<+% taxonomy.tip +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.1))) +
  geom_treescale(x = 0, y = 0) +
  geom_tiplab(size = 4, hjust = -0.1, aes(label = taxon)) +
  geom_nodelab(aes(label = label), size = 2, hjust = -0.3); mix.tip.tree.collapsed
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/Mitochondrial_phylogeny.pdf", plot = mix.tip.tree.collapsed, width = 12, height = 16, dpi = 1200)
