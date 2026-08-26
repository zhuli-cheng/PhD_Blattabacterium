library(Biostrings)

source("Main.R")

#read full fasta
SSU.all <- readDNAStringSet("../input/16S.txt", format="fasta")
names(SSU.all)

#decide samples to exclude
families.to.PCR <- c("Blattellidae", "Pseudophyllodromiidae", "Anaplectidae")
families.all <- c("Blattellidae", "Blaberidae", "Pseudophyllodromiidae", "Ectobiidae", "Anaplectidae", "Blattidae", "Corydiidae")
families.to.exclude <- setdiff(families.all, families.to.PCR)
samples.to.exclude <- taxonomy %>%
  filter(family %in% families.to.exclude | family_inferred %in% families.to.exclude)
to.exclude <- paste(samples.to.exclude$sample, collapse = "|") 
to.exclude <- paste(to.exclude, "|ISG-1_S63_L003|YS-1fin_S58_L003|OKE030_S43_L003|AU023_S352") #removing some sequences that cause gaps

#exclude and save
SSU.filtered <- SSU.all[!grepl(to.exclude, names(SSU.all))]
writeXStringSet(SSU.filtered, filepath = "../output/16S_filtered.fasta", format = "fasta")
names(SSU.filtered)
