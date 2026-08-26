library(dplyr)
library(ggplot2)

######### length table for kinjo-san #######
length.table <- dataset.clean %>%
  group_by(orthogroups_in_OG) %>%
  mutate(gene=paste(toString(unique(Preferred_name)), toString(unique(Description)))) %>%
  ungroup() %>%
  select(orthogroups_in_OG, query, position, nucleotides, GC, type, curated_COG_category, OGs_in_orthogroup, gene)

length(unique(length.table$orthogroups_in_OG))
length(unique(length.table$OGs_in_orthogroup))

#check orthogroups w/o COG categories
COG.empty <- length.table %>%
  group_by(orthogroups_in_OG) %>%
  mutate(COG_categories = toString(unique(na.omit(curated_COG_category)))) %>%
  ungroup() %>%
  filter(COG_categories == "")
COG.empty <- as.character(COG.empty$OGs_in_orthogroup)
View(dataset[dataset$OGs_in_orthogroup %in% COG.empty, ])

#all with unknown functions, so give them category S
condition <- length.table$OGs_in_orthogroup %in% COG.empty
length.table$curated_COG_category[condition] <- "S"

#plot number of orthogroups in each COG category
ggplot(length.table, aes(x = curated_COG_category)) +
  geom_bar()

#summarize the average length in each orthogroup
length.table.summary <- length.table %>%
  group_by(orthogroups_in_OG) %>%
  mutate(mean_length = mean(na.omit(nucleotides))) %>%
  mutate(COG_categories = names(which.max(table(curated_COG_category)))) %>%
  ungroup() %>%
  select(orthogroups_in_OG, mean_length, COG_categories, OGs_in_orthogroup, gene) %>%
  unique()

#write.csv(length.table, "length.table.csv")
#write.csv(length.table.summary, "length.table.summary.csv")

