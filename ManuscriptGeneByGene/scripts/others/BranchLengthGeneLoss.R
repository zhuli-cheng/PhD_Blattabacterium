#source("Main.R")
#source("CollapsingHighlySimilarGenomes.R")
source("AncestralLossReconstruction.R")

#calculate root-to-tip distance
depths <- node.depth.edgelength(collapsed.tree)
root.to.tip <- depths[1:Ntip(collapsed.tree)]
names(root.to.tip) <- collapsed.tree$tip.label

#correlation between root-to-tip and gene loss
distance.loss <- sample.orthogroup %>%
  select(sample, mean_GC_of_orthologs, no_of_orthologs, family)

distance.loss$root_to_tip <- root.to.tip[distance.loss$sample]

p2 <- ggplot(data = distance.loss, aes(y = max(no_of_orthologs) - no_of_orthologs, x = root_to_tip, color = family)) + 
  geom_point() + 
  xlab(expression("Root-to-tip distance of a "~italic("Blattabacterium")~" genome on the phylogeny")) +
  ylab("Number of lost PCGs in the genome") +
  labs(color = "Host") +
  theme_classic(); p2

ggsave("~/Desktop/Thesis/Typst/figures/chap3results/Distance_gene_loss.png", plot = p2, width = 10, height = 5, dpi = 1500)

#heatmap
rownames(distance.loss) <- distance.loss$sample

distance.loss <- distance.loss %>%
  mutate(gene_loss = max(no_of_orthologs) - no_of_orthologs) %>%
  mutate(GC_loss = max(mean_GC_of_orthologs) - mean_GC_of_orthologs) %>%
  select(root_to_tip, GC_loss, gene_loss)

png("~/Desktop/Thesis/Typst/figures/chap3results/gene_loss_heatmap.png", width = 1600, height = 1000)  # Adjust the size as needed
phylo.heatmap(collapsed.tree, distance.loss, fsize = c(0.01, 1.8) , standardize = TRUE, 
              legend = TRUE, labels = TRUE, split = c(0.8, 0.2))
dev.off()



