
gc.heatmap <- function(orthogroups, pathway.description) {
  #gene description text
  gene.description <- orthogroup.sample %>%
    filter(new_orthogroup %in% orthogroups) %>%
    mutate(new_orthogroup = factor(new_orthogroup, levels = orthogroups)) %>%
    arrange(new_orthogroup) %>%
    transmute(`Genes in the pathway` = str_wrap(gene, width = 140))
  
  #gene GC content dataset
  subset <- sample.orthogroup %>%
    select(paste0("new_orthogroup", orthogroups))
  rownames(subset) <- sample.orthogroup$sample
  subset[]<- lapply(subset, as.numeric)
  subset <- subset[match(collapsed.tree$tip.label, rownames(subset)), ]
  subset[subset== 0] <- NA
  
  #plot
  png("heatmap_plot.png", width = 1600, height = 1000)  # Adjust the size as needed
  phylo.heatmap(collapsed.tree, subset, fsize = 0.001, standardize = FALSE, 
                legend = TRUE, labels = TRUE, split = c(0.8, 0.2))
  title(pathway.description, line = -3, cex.main = 4)
  dev.off()  # Close the device to save the image
  heatmap_img <- rasterGrob(readPNG("heatmap_plot.png"), width = unit(1, "npc"), height = unit(1, "npc"))
  
  # Create the table grob
  table_grob <- tableGrob(
    gene.description,
    rows = seq_len(nrow(gene.description)),
    theme = ttheme_default(
      base_size = 16,
      core = list(fg_params = list(hjust = 0, x = 0)),
      colhead = list(fg_params = list(hjust = 0, x = 0))
    )
  )
  table_grob$widths <- unit(c(0.05, 0.90), "npc")
  
  # Combine the heatmap image and table
  grid.arrange(heatmap_img, table_grob, nrow = 2, heights = c(0.7, 0.3))
}



