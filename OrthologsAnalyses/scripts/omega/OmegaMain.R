source("Main.R")
source("pathways/KEGGtable.R")

########### omega ~ number of families with gene loss
M0.omega <- read.csv(paste0("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsAnalyses/", sample.size, "/input/omega_summary.out"), header = T)
family.omega <- data.frame(tip_label = collapsed.tree$tip.label) %>%
  mutate(host_family = taxonomy$family[match(collapsed.tree$tip.label, taxonomy$sample)]) 

#create table showing the loss of genes in each families
combined_family_presence <- data.frame()
for (i in unique(dataset.clean$new_orthogroup)) {
  #load gene presence/absence
  orthogroup.name <- paste0("new_orthogroup", i)
  orthogroup.presence <- setNames(ifelse(sample.orthogroup[[orthogroup.name]] == 0, 0, 1), sample.orthogroup$sample)
  
  #find families that have gene loss
  family.presence <- data.frame(tip_label = collapsed.tree$tip.label) %>%
    mutate(family = taxonomy$family[match(collapsed.tree$tip.label, taxonomy$sample)]) %>%
    mutate(gene_presence = orthogroup.presence[collapsed.tree$tip.label]) %>%
    group_by(family) %>%
    mutate(family_loss = ifelse(sum(gene_presence) == 0, "complete_gene_loss", ifelse(sum(gene_presence) == n(), "no_gene_loss", "partial_gene_loss"))) %>%
    mutate(family_with_loss = ifelse(sum(gene_presence) == n(), 0, 1)) %>%
    ungroup() %>%
    select(family, family_loss, family_with_loss) %>%
    unique() %>%
    mutate(number_family_with_loss = sum(family_with_loss)) %>%
    mutate(new_orthogroup = i) 
  
  combined_family_presence <- rbind(combined_family_presence, family.presence)
}

combined_family_presence <- combined_family_presence %>%
  select(new_orthogroup, number_family_with_loss) %>%
  unique()

M0.omega <- full_join(combined_family_presence, M0.omega, by = c("new_orthogroup")) %>%
  full_join(orthogroup.sample, by = "new_orthogroup") %>%
  select(new_orthogroup, number_family_with_loss, omega, no_of_samples, curated_COG_category, gene, n, dN, dN_standardized, dS, dS_standardized)

#M0.omega <- M0.omega %>%
#  filter(omega < 2) 

M0.omega <- left_join(M0.omega, orthogroup.sample[c("new_orthogroup", "mean_position_in_orthogroup", "mean_GC_in_orthogroup", "sd_GC_in_orthogroup")], by = c("new_orthogroup"))

#combine KEGG information
ko.categories <- ko.table %>%
  select(!Preferred_name) %>%
  separate_rows(new_orthogroup, sep = ", ") %>%
  mutate(new_orthogroup = as.numeric(new_orthogroup)) %>%
  left_join(M0.omega, by = "new_orthogroup")

ko.categories <- ko.categories %>%
  mutate(original_order = row_number()) %>%  # Store original order
  filter(pathway %in% pathway.list) %>%
  arrange(original_order) %>%  # Restore original order
  mutate(category2 = factor(category2, levels = unique(category2))) %>%
  as.data.frame()


#split double COG categories into single duplicated ones
M0.omega.COG <- M0.omega %>%
  mutate(curated_COG_category = str_replace_all(curated_COG_category, ",", "")) %>%
  mutate(curated_COG_category = str_replace_all(curated_COG_category, " ", "")) %>%
  mutate(curated_COG_category = str_split(curated_COG_category, "")) %>%
  unnest(curated_COG_category)


