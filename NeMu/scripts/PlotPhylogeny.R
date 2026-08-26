source("Main.R")

subject.mapping.list <- list()

for (group in selected.groups$group) {
  tree <- read.tree(paste0("../input/monophyleticGroups/NeMuPipeline/", group, "/output/final_tree.nwk"))
  ggtree(tree) +
    geom_tiplab(size = 2.5, hjust = -0.1) +  
    geom_treescale(x = 0, y = 0) +
    scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +
    geom_nodelab()
  #load subject mapping file
  subject.mapping <- read.table(paste0("../input/monophyleticGroups/NeMuPipeline/", group, "/output/species_mapping.txt"), header = FALSE, sep = "")
  subject.mapping.list[[group]] <- subject.mapping
}


group <- "group2"


