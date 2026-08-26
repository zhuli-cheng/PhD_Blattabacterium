#source("Main.R")

dataset.kegg <- dataset.clean %>%
  select(new_orthogroup, Preferred_name, Description, EC, KEGG_ko, KEGG_Pathway, mean_position_in_orthogroup) %>%
  unique() %>%
  group_by(new_orthogroup) %>%
  summarise(across(c(Preferred_name, Description, EC, KEGG_ko, KEGG_Pathway), ~ paste(., collapse = ", ")),
            mean_position_in_orthogroup = first(mean_position_in_orthogroup)) %>%
  ungroup() 

dataset.kegg$KEGG_Pathway[dataset.kegg$Preferred_name == "sad"] <- "NA"

dataset.kegg <- dataset.kegg[order(dataset.kegg$mean_position_in_orthogroup), ]


ko.table <- read.csv("../input/KEGG_table.csv", header = F, colClasses = c("V3" = "character"))
ko.table$V3 <- paste0("ko", ko.table$V3)
colnames(ko.table) <- c("category1", "category2", "ko", "pathway")

kaas.table <- read.table("../input/KAAS.ko.txt", na.strings="", fill=TRUE)
colnames(kaas.table) <- c("new_orthogroup", "kaas_ko")

#write.table(
#  kaas.table,
#  "~/Desktop/Thesis/Typst/figures/supplementary materials//KO_assignments.csv",
#  sep = ",",
#  row.names = FALSE,
#  col.names = TRUE,
#  quote = TRUE)


tmp <- full_join(dataset.kegg, kaas.table, by = join_by(new_orthogroup))

find_genes <- function(ko) {
  match_row <- dataset.kegg %>%
    filter(grepl(ko, KEGG_Pathway)) %>% # Check if KO exists in ko_numbers column
    pull(Preferred_name)
  
  if (length(match_row) == 0) return(NA) # Return NA if no match found
  return(match_row)
}
ko.table$Preferred_name <- sapply(ko.table$ko, function(x) paste(find_genes(x), collapse = ", "))

find_orthogroups <- function(ko) {
  match_row <- dataset.kegg %>%
    filter(grepl(ko, KEGG_Pathway)) %>% # Check if KO exists in ko_numbers column
    pull(new_orthogroup)
  
  if (length(match_row) == 0) return(NA) # Return NA if no match found
  return(match_row)
}
ko.table$new_orthogroup <- sapply(ko.table$ko, function(x) paste(find_orthogroups(x), collapse = ", "))

ko.table <- ko.table %>%
  filter(new_orthogroup != "NA")
  
#write.table(
#  ko.table,
#  "~/Desktop/Thesis/Typst/figures/supplementary materials//KEGG_pathways.csv",
#  sep = ",",
#  row.names = FALSE,
#  col.names = TRUE,
#  quote = TRUE)


tmp <- ko.table %>%
  select(category1, category2, pathway, Preferred_name)

#write.csv(tmp, "../244/output/ko_table.csv", row.names = F)

pathway.list <- c(
  "Glycolysis_or_Gluconeogenesis",
  "Citrate_cycle_TCA_cycle",
  "Pentose_phosphate_pathway",
  "Pyruvate_metabolism",
  "Glyoxylate_and_dicarboxylate_metabolism",
  "C5_Branched_dibasic_acid_metabolism",
  "Oxidative_phosphorylation",
  "Nitrogen_metabolism",
  "Sulfur_metabolism",
  "Fatty_acid_biosynthesis",
  "Glycerolipid_metabolism",
  "Glycerophospholipid_metabolism",
  "Purine_metabolism",
  "Pyrimidine_metabolism",
  "Alanine_aspartate_and_glutamate_metabolism",
  "Glycine_serine_and_threonine_metabolism",
  "Cysteine_and_methionine_metabolism",
  "Valine_leucine_and_isoleucine_biosynthesis",
  "Lysine_biosynthesis",
  "Arginine_biosynthesis",
  "Histidine_metabolism",
  "Phenylalanine_tyrosine_and_tryptophan_biosynthesis",
  "Peptidoglycan_biosynthesis",
  "Thiamine_metabolism",
  "Riboflavin_metabolism",
  "Nicotinate_and_nicotinamide_metabolism",
  "Pantothenate_and_CoA_biosynthesis",
  "Vitamin_B6_metabolism",
  "Biotin_metabolism",
  "Folate_biosynthesis",
  "One_carbon_pool_by_folate",
  "Lipoic_acid_metabolism",
  "Porphyrin_metabolism",
  "Ubiquinone_and_other_terpenoid_quinone_biosynthesis",
  "Terpenoid_backbone_biosynthesis",
  "RNA_polymerase",
  "Ribosome",
  "Aminoacyl_tRNA_biosynthesis",
  "Protein_export",
  "RNA_degradation",
  "DNA_replication",
  "Base_excision_repair",
  "Mismatch_repair",
  "Homologous_recombination",
  "ABC_transporters",
  "Bacterial_secretion_system",
  "Two_component_system"
)



#write.table(
#  gsub("_", " ", pathway.list), 
#  "~/Desktop/Thesis/Typst/figures/chap2methods/KEGG_selection.csv",
#  sep = ",",
#  row.names = FALSE,
#  col.names = FALSE,
#  quote = TRUE)


