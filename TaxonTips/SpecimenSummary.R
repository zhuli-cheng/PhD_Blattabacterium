library(dplyr)
library(readxl)
library(tidyr)
library(gt)

specimen <- read_excel("~/Desktop/Zhuli_BLT446.xlsx", sheet = "main") %>%
  filter(Sample != "1N-T11_S95_L002") %>%
  filter(Sample != "1N-Z1909_S96_L002")

taxonomy <- read.csv("~/Desktop/BlattabacteriumRScripts/TaxonTips/input/host.csv")
colnames(taxonomy) <- c("Sample", "family", "family_inferred", "species")
taxonomy <- taxonomy %>%
  mutate(family = ifelse(family == "", family_inferred, family)) %>%
  mutate(mixed = ifelse(species == "", family, species)) 
  
specimen.taxonomy <- left_join(specimen, taxonomy, by = c("Sample")) %>%
  mutate(Species = mixed)

specimen.appendix <- specimen.taxonomy %>%
  select(Sample, Species, "Extraction date", Country, Location, GPS)

write.csv(file="~/Desktop/Thesis/Typst/figures/appendices/specimens446.csv", specimen.appendix, row.names = FALSE)

list.106 <- read.table("input/list106.txt")
specimen.manuscript <- specimen.appendix %>%
  filter(Sample %in% list.106$V1)

write.csv(file="~/Desktop/Thesis/Typst/figures/supplementary materials/specimens106.csv", specimen.manuscript, row.names = FALSE)


