#source("Main.R")
#source("CollapsingHighlySimilarGenomes.R")
source("AncestralLossReconstruction.R")

################# show each gene's GC content of a gene on the tree ###############
#specify orthogroups to look at
GC_on_tree <- function(choose_orthogroup) {
  orthogroup_name <- paste0("orthogroup", choose_orthogroup)
  assign(orthogroup_name, setNames(as.numeric(sample.orthogroup[[paste0("new_orthogroup", choose_orthogroup)]]), sample.orthogroup$sample))
  
  # Check if all names are present in the tree
  orthogroup_data <- get(orthogroup_name)
  if (!all(names(orthogroup_data) %in% tree$tip.label) || !all(tree$tip.label %in% names(orthogroup_data))) {
    stop("Mismatch between the tree tip labels and the orthogroup data names.")
  }
  
  collapsed.tree$tip.label <- paste0(collapsed.tree$tip.label, ":", orthogroup_data[collapsed.tree$tip.label])
  
  tip.data <- data.frame(tip_label = collapsed.tree$tip.label) %>%
    mutate(condition = ifelse(grepl(":0", tip_label), "gene lost", "gene retained"))
  
  orthogroup_GC_plot <- ggtree(collapsed.tree) %<+% tip.data +
    geom_tiplab(size = 2.5, hjust = -0.1, aes(color = condition)) +
    scale_color_manual(values = c("gene lost" = "red", "gene retained" = "black")) +
    theme(plot.margin = margin(0, 0, 0, 0, "cm"), legend.position.inside = c(0, 1), 
          legend.justification = c("left", "top"), legend.title = element_blank())
  
  print(orthogroup_GC_plot)
  # ggsave(filename = paste0("orthogroup", choose_orthogroup, "dnaN.GC.png"), plot = orthogroup_GC_plot, width = 10, height = 15, dpi = 1200)
}

GC_on_tree(1)


################### plot gene loss ~ GC in genomes ################
p1 <- ggplot(data = sample.orthogroup, aes(x = mean_GC_of_orthologs, y = max(no_of_orthologs) - no_of_orthologs, color = family)) + 
  geom_point() + 
  xlab(expression("average GC-content of PCGs in a "~italic("Blattabacterium")~" genome")) +
  ylab(expression("Number of lost PCGs in the genome")) +
  labs(color = "Host") +
  theme_classic(); p1

#ggsave("~/Desktop/Thesis/Typst/figures/chap3results/Genome_GC_gene_loss.png", plot = p1, width = 10, height = 5, dpi = 1500)

################### plot gene loss ~ GC in orthogroups ################
GC.gene.loss <- gene.loss %>%
  full_join(orthogroup.sample, by = "new_orthogroup")

plot(data = GC.gene.loss, gene_loss ~ mean_GC_in_orthogroup)
summary(lm(data = GC.gene.loss, gene_loss ~ mean_GC_in_orthogroup))

plot(data = GC.gene.loss, gene_loss ~ sd_GC_in_orthogroup)
summary(lm(data = GC.gene.loss, gene_loss ~ sd_GC_in_orthogroup))

plot(data = GC.gene.loss, gene_loss ~ mean_length_in_orthogroup)
summary(lm(data = GC.gene.loss, gene_loss ~ mean_length_in_orthogroup))


p <- ggplot(data = GC.gene.loss, aes( x = mean_GC_in_orthogroup, y = gene_loss, color = curated_COG_category)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method = "lm", se = FALSE, aes(group = curated_COG_category)) 

ggplot(data = GC.gene.loss, aes( x = sd_GC_in_orthogroup, y = gene_loss, color = curated_COG_category)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = curated_COG_category))   # Add linear model

ggplot(data = GC.gene.loss, aes( x = mean_length_in_orthogroup, y = gene_loss, color = curated_COG_category)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, aes(group = curated_COG_category))   # Add linear model

