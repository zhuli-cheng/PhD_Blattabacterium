library(ggplot2)
library(reshape2)
library(lme4)
library(dplyr)
library(phytools)

gene_GC_genome_GC <- function(choose_orthogroup) {
  data <- dataset.clean %>%
    subset(new_orthogroup == choose_orthogroup) %>%
    select(query, sample, nucleotides, GC, type) %>%
    right_join(sample.summary, by = "sample") %>%
    mutate(GC = ifelse(is.na(GC), 0, GC))
  
  plot <- ggplot(data, aes(x = mean_GC_of_orthologs, y = GC, color = type)) +
    geom_point() +
    ylim(0, 40) +
    theme_minimal() +
    labs(x = "genome GC", y = "gene GC", title = paste0("new orthogroup ", choose_orthogroup)) 
  
  print(plot)
}

for (group in 599:500) {
  gene_GC_genome_GC(group)
}


gene_GC_genome_GC(505)


#linear models between GC contents
summary(lm(new_orthogroup507 ~ mean_GC_of_orthologs, data = sample.orthogroup))

summary(lm(new_orthogroup..1 ~ new_orthogroup.13, data = sample.orthogroup))

################## plot gene gc to genome gc heatmap ############################
#specify orthogroups to look at
gene597 <- sample.orthogroup %>%
  select(no_of_orthologs, mean_GC_of_orthologs, new_orthogroup597)
gene597[]<- lapply(gene597, as.numeric)
rownames(gene597) <- sample.orthogroup$sample

gene597 <- gene597[match(collapsed.tree$tip.label, rownames(gene597)), ]

gene597[gene597== 0] <- NA

phylo.heatmap(collapsed.tree, gene597, fsize = 0.001, standardize = TRUE, legend = FALSE, labels = FALSE, split=c(0.8, 0.2))
#add.scale.bar(cex = 0.1, pos = "bottomleft")

#larger one to crop the unit bar
phylo.heatmap(collapsed.tree, gene597, fsize = 1, standardize = TRUE, legend = TRUE, labels = FALSE, split=c(2, 8))



