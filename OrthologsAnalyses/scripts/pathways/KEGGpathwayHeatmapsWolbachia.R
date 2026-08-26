#load mito tree 
setwd("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/MitochondrialTree/scripts/")
source("OtherSymbionts.R")

#set working directory back
setwd("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsAnalyses/scripts/")
source("pathways/KEGGtable.R")


collapsed.tree <- mito.tree
#no_blattabacterium <- rownames(symbionts[symbionts$absence_confirmed == 1, ])
#collapsed.tree <- drop.tip(collapsed.tree, no_blattabacterium)
#collapsed.tree <- drop.tip(collapsed.tree, c("CAM19_2_C164_CR_S63_L004", "CAM19_1_C037_CR_S37_L001", "CAM19_2_C007_CR_S9_L002"))



wolbachia.heatmap <- function(pathway.description) {
  orthogroups <- ko.table %>%
    filter(pathway == pathway.description | category2 == pathway.description | category1 == pathway.description) %>%
    pull(new_orthogroup) %>%
    strsplit(", ") %>%
    unlist() %>%
    unique()
  
  #gene GC content dataset
  subset <- sample.orthogroup %>%
    select(paste0("new_orthogroup", orthogroups))
  rownames(subset) <- sample.orthogroup$sample
  subset[]<- lapply(subset, as.numeric)
  subset <- subset[rownames(subset) %in% collapsed.tree$tip.label, ]
  
  collapsed.tree <- keep.tip(collapsed.tree, rownames(subset))
  
  subset <- subset[match(collapsed.tree$tip.label, rownames(subset)), ]
  subset[subset == 0] <- NA
  subset[subset > 0] <- 50
  
  
  #add Wolbachia and Rickettsia
  symbionts$wolbachia_presence[symbionts$wolbachia_presence > 0] <- 100
  symbionts$rickettsia_presence[symbionts$rickettsia_presence > 0] <- 70
  
  subset <- subset %>%
    mutate(Wolbachia = symbionts$wolbachia_presence[match(rownames(subset), rownames(symbionts))]) %>%
    mutate(Rickettsia = symbionts$rickettsia_presence[match(rownames(subset), rownames(symbionts))])
  
  
  #plot
  png("heatmap_plot.png", width = 1600, height = 1000)  # Adjust the size as needed
  phylo.heatmap(collapsed.tree, subset, fsize = 0.001, standardize = FALSE, 
                legend = FALSE, labels = TRUE, split = c(0.8, 0.2))
  title(pathway.description, line = -3, cex.main = 4)
  dev.off()  # Close the device to save the image
  heatmap_img <- rasterGrob(readPNG("heatmap_plot.png"), width = unit(1, "npc"), height = unit(1, "npc"))
  
  # legends
  gene.description <- orthogroup.sample %>%
    filter(new_orthogroup %in% orthogroups) %>%
    mutate(new_orthogroup = factor(new_orthogroup, levels = orthogroups)) %>%
    arrange(new_orthogroup) %>%
    transmute("Genes and Symbionts" = str_wrap(gene, width = 140))
  
  gene.description[nrow(gene.description) + 1, ] <- "Wolbachia contigs"
  gene.description[nrow(gene.description) + 1, ] <- "Rickettsia contigs"

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

kegg.pathway <- function(pathway.description) {
  png(paste0("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsAnalyses/", sample.size, "/output/KEGG_pathways_Wolbachia/", pathway.description, ".png") , width = 1200, height = 1600)
  wolbachia.heatmap(pathway.description)
  dev.off() 
}

#Metabolism of cofactors and vitamins
#kegg.pathway("Metabolism of cofactors and vitamins") #category2
kegg.pathway("Riboflavin_metabolism") #B2
kegg.pathway("Vitamin_B6_metabolism") #B6
kegg.pathway("Biotin_metabolism") #B7
kegg.pathway("Folate_biosynthesis") #B9

