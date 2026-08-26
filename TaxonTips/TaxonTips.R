library(dplyr)
library(readxl)

#modify outgroup tip
outgroups <- read.csv("~/Desktop/RScripts/TaxonTips/input/outgroup_species.csv", header = TRUE)
colnames(outgroups) <- c("sample", "species") 
outgroups <- outgroups %>%
  mutate(species = gsub(" ", "_", species)) %>%
  mutate(family = species) %>%
  mutate(mixed = species)

#process taxonomy file
taxonomy <- read.csv("~/Desktop/RScripts/TaxonTips/input/host.csv")
colnames(taxonomy) <- c("sample", "family", "family_inferred", "species")
taxonomy <- taxonomy %>%
  mutate(family = ifelse(family == "", family_inferred, family)) %>%
  mutate(mixed = ifelse(species == "", family, species)) 

taxonomy <- bind_rows(taxonomy, outgroups)
rm(list = "outgroups")

#use functions to make new files to replace tip_labels
#can use "family", "species", "mixed" as level
modify_tip <- function(tree, level){
  taxonomy <- taxonomy %>%
    select(sample, taxon = all_of(level))
  taxonomy.tip <<- data.frame(tip_label = tree$tip.label) %>%
    mutate(taxon = sapply(tree$tip.label, function(tip) {
      match_idx <- which(sapply(taxonomy$sample, function(sample) grepl(sample, tip)))
      taxonomy$taxon[match_idx] })) %>%
      mutate(taxon = ifelse(taxon == "", tip_label, taxon)) %>%
      mutate(taxon = sapply(tree$tip.label, function(tip) {
        match_idx <- which(sapply(taxonomy$sample, \(s) grepl(s, tip)))
        if (length(match_idx)) taxonomy$taxon[match_idx[1]] else tip}))
}


