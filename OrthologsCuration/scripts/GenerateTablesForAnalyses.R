library(dplyr)
library(tidyr)
library(tibble)

#list of all orthologs
curated.orthologs <- dataset.clean %>%
  select(query, new_orthogroup, curated_COG_category, OGs_in_orthogroup, sample_size_in_orthogroup, type, sample, GC, nucleotides)

#make sure there is no , in the elements, then
#write.csv(curated.orthologs, paste0("../", sample.size, "/output/curated_orthologs.csv"))

#sample ~ orthogroup table
orthogroup <- dataset.clean %>%
  group_by(orthogroups_in_OG, OGs_in_orthogroup, new_orthogroup) %>%
  mutate(gene=paste(toString(unique(Preferred_name)), toString(unique(Description)))) %>%
  ungroup() %>%
  select(new_orthogroup, curated_COG_category, OGs_in_orthogroup, gene, mean_position_in_orthogroup, mean_GC_in_orthogroup, sd_GC_in_orthogroup, mean_length_in_orthogroup) %>%
  count(new_orthogroup, curated_COG_category, OGs_in_orthogroup, gene, mean_position_in_orthogroup, mean_GC_in_orthogroup, sd_GC_in_orthogroup, mean_length_in_orthogroup)

orthogroup.sample <- dataset.clean %>%
  select(new_orthogroup, sample, GC) %>%
  pivot_wider(names_from = sample, values_from = GC) %>%
  full_join(orthogroup, by = join_by(new_orthogroup)) %>%
  select(new_orthogroup, curated_COG_category, OGs_in_orthogroup, gene, n, everything()) %>%
  as.data.frame


write.csv(orthogroup.sample, paste0("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsCuration/", sample.size, "/output/orthogroup_sample.csv"), row.names = F)


# transpose the dataframe to use in ape
sample.orthogroup <- t(orthogroup.sample)
colnames(sample.orthogroup) <- paste0("new_orthogroup", as.numeric(sample.orthogroup["new_orthogroup", ]))

#check if there is "," in sample.orthogroup. If not, run: 
sample.orthogroup <- as.data.frame(sample.orthogroup)
sample.orthogroup <- sample.orthogroup[grepl("_S", rownames(sample.orthogroup)), ]
sample.orthogroup[] <- lapply(sample.orthogroup, as.numeric)
sample.orthogroup$mean_GC_of_orthologs <-  rowMeans(sample.orthogroup, na.rm = TRUE)
sample.orthogroup$no_of_orthologs <- rowSums(!is.na(sample.orthogroup)) - 1
sample.orthogroup <- sample.orthogroup %>%
  rownames_to_column(var = "sample")
  
write.csv(sample.orthogroup, paste0("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsCuration/", sample.size, "/output/sample_orthogroup.csv"), row.names = F)


write.csv(dataset.clean, paste0("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsCuration/", sample.size, "/output/dataset_clean.csv"), row.names = F)

#table containing headers 
orthogroup.headers <- dataset.clean %>%
  select(query, new_orthogroup)
write.table(orthogroup.headers, paste0("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsCuration/", sample.size, "/output/orthogroup_headers.csv"), row.names = F, col.names = F, quote = F, sep = ",")


