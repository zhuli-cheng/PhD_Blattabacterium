######### GC heatmap
GC <- GC.absence %>%
  mutate(blattabacterium = as.numeric(25)) 
GC <- GC[mito.tree$tip.label, , drop = FALSE]
phylo.heatmap(mito.tree, GC[, c("blattabacterium", "GC")], fsize = 0.01, standardize = TRUE, legend = FALSE, labels = FALSE, split=c(0.8, 0.2)) 

############### number of orthologs
#specify orthogroups to look at
no.of.orthogroups <- setNames(as.numeric(sample.orthogroup$no_of_orthologs), sample.orthogroup$sample)
all(names(no.of.orthogroups) %in% mito.tree$tip.label)
all(mito.tree$tip.label %in% names(no.of.orthogroups))

no.of.orthogroups <- no.of.orthogroups[rownames(symbionts)]
symbionts$gene.repertoire <- no.of.orthogroups

#map statistics on the tree tips
gene.repertoire.tree <- mito.tree
gene.repertoire.tree$tip.label <- paste0(gene.repertoire.tree$tip.label, ":", symbionts$gene.repertoire[match(gene.repertoire.tree$tip.label, rownames(symbionts))])

tip.data <- data.frame(tip_label = gene.repertoire.tree$tip.label,
                       blattabacterium_absence = absence$blattabacterium_absence[match(mito.tree$tip.label, absence$sample)])
rownames(tip.data) <- mito.tree$tip.label

gene.repertoire.loss <- ggtree(gene.repertoire.tree) %<+% tip.data +
  geom_tiplab(size = 3, hjust = 0, aes(color = blattabacterium_absence)) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), legend.position = c(0, 1), legend.justification = c("left", "top"), legend.title = element_blank())
ggsave(filename = "../output/gene_repertoire_loss.pdf", plot = gene.repertoire.loss, width = 30, height = 12, dpi = 3600)

#remove incomplete genomes, to only keep circular complete genomes to show individual gene losses
tip.data.circular <- tip.data %>%
  filter(!grepl(":NA", tip_label) | blattabacterium_absence == 1) 
tip.to.drop <- setdiff(gene.repertoire.tree$tip.label, tip.data.circular$tip_label)
gene.repertoire.tree.circular <- drop.tip(gene.repertoire.tree, tip.to.drop)
gene.repertoire.loss.circular <- ggtree(gene.repertoire.tree.circular) %<+% tip.data.circular +
  geom_tiplab(size = 1, hjust = -0.1, aes(color = blattabacterium_absence)) +
  theme(plot.margin = margin(1, 1, 1, 1, "cm"), legend.position = c(0, 1), legend.justification = c("left", "top"), legend.title = element_blank())
ggsave(filename = "../output/gene_repertoire_loss_circular.pdf", plot = gene.repertoire.loss.circular, width = 30, height = 12, dpi = 3600)

