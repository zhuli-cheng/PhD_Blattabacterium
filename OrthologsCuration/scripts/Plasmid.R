#genes that are transferred between chromosome/plasmid
chromosome.plasmid <- dataset %>%
  subset (locations_in_orthogroup == "plasmid, chromosome")

unique(chromosome.plasmid$Description)
unique(chromosome.plasmid$orthogroup)

plot((chromosome.plasmid$position[chromosome.plasmid$location == "chromosome"]))
#these genes are found on three locations in the chromosome

#integration of plasmid onto the chromosome? 
plasmid <- dataset %>%
  subset(location == "plasmid") 

unique(plasmid$Description)
unique(plasmid$orthogroup)

plasmid.genes <- plasmid %>%
  select(Preferred_name, EC, KEGG_ko, PFAMs, Description) %>%
  unique() %>%
  filter(!is.na(Description))

write.table(
  plasmid.genes,
  "~/Desktop/Thesis/Typst/figures/chap3results/plasmid_genes.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = TRUE)
