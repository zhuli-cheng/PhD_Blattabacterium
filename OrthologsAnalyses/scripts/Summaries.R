source("LoadPackages.R")
source("ChromosomePosition.R")

orthogroups <- orthogroup.sample %>%
  select(new_orthogroup, n, curated_COG_category, OGs_in_orthogroup, gene) 

#order the orthogroups by their relative position on the chromosome
orthogroups_ordered <- orthogroups[match(position$new_orthogroup, orthogroups$new_orthogroup), ]

orthogroups_ordered <- orthogroups_ordered %>%
  mutate(position = row_number()) %>%
  select(position, n, curated_COG_category, OGs_in_orthogroup, gene)

write.table(
  orthogroups_ordered,
  "~/Desktop/Thesis/Typst/figures/supplementary materials/orthogroups.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = TRUE)
