source("LoadPackages.R")
source("ThemeAnimalBehaviour.R")

source("../../TaxonTips/TaxonTips.R")

####################### list of samples #################################
list <- read.table("../input/list")

###################### full data on each protein-coding gene #############
dataset.clean <- read.csv("/Users/zhulicheng/Desktop/RScripts/OrthologsCuration/244/output/dataset_clean.csv", header = T) %>%
  filter(sample %in% list$V1) %>%
  group_by(new_orthogroup) %>%
  mutate(genes_in_orthogroup = n()) %>%
  mutate(sample_size_in_orthogroup = n_distinct(sample)) %>%
  mutate(mean_position_in_orthogroup = mean(na.omit(position))) %>%
  mutate(mean_GC_in_orthogroup = mean(na.omit(GC))) %>%
  mutate(sd_GC_in_orthogroup = sd(na.omit(GC))) %>%
  mutate(mean_length_in_orthogroup = mean(na.omit(nucleotides))) %>%
  mutate(median_length_in_orthogroup = median(na.omit(nucleotides))) %>%
  mutate(relative_length_in_orthogroup = nucleotides/median_length_in_orthogroup) %>%
  mutate(distance_length_in_orthogroup = nucleotides - median_length_in_orthogroup) %>%
  mutate(max_length_in_orthogroup = max(na.omit(relative_length_in_orthogroup))) %>%
  mutate(min_length_in_orthogroup = min(na.omit(relative_length_in_orthogroup))) %>%
  mutate(OGs_in_orthogroup = ifelse(length(unique(curated_OG[!is.na(curated_OG)])) == 0, NA, toString(unique(curated_OG[!is.na(curated_OG)])))) %>% 
  mutate(OGs_in_orthogroup = ifelse(!is.na(OGs_in_orthogroup), OGs_in_orthogroup, ifelse(!is.na(unique(OG_root[!is.na(OG_root)])), toString(unique(OG_root[!is.na(OG_root)])), NA))) %>%
  mutate(locations_in_orthogroup = toString(unique(location))) %>%
  ungroup() %>%
  mutate(OG_inferred_from_orthogroup = ifelse(!is.na(curated_OG), curated_OG, OGs_in_orthogroup)) %>%
  group_by(OG_inferred_from_orthogroup) %>%
  mutate(orthogroups_in_OG = ifelse(is.na(OG_inferred_from_orthogroup), NA, toString(unique(new_orthogroup)))) %>%
  ungroup() %>%
  mutate(core_proteome = ifelse(sample_size_in_orthogroup == 81, TRUE, FALSE)) %>%
  group_by(new_orthogroup) %>%
  mutate(curated_COG_category = toString(unique(na.omit(curated_COG_category)))) %>%
  ungroup() 

################ by sample ######################
sample.orthogroup <- read.csv("/Users/zhulicheng/Desktop/RScripts/OrthologsCuration/244/output/sample_orthogroup.csv", header = T) %>%
  filter(sample %in% list$V1)

sample.orthogroup[is.na(sample.orthogroup)] <- "0"
sample.orthogroup$family <- taxonomy$family[match(sample.orthogroup$sample, taxonomy$sample)]
sample.summary <- sample.orthogroup %>%
  select(sample, no_of_orthologs, mean_GC_of_orthologs)

################# tree #########################
tree <- read.tree("../input/81_RY.treefile")

#reroot the tree
two.families <- c("CAM19_2_C149_CR_S49_L004", "YS-3-1_S60_L003")
common.ancestor <- getMRCA(tree, two.families); common.ancestor
tree <- root(tree, node = common.ancestor, resolve.root = T)
tree$node.label[tree$node.label == "Root"] <- ""

#visualize the tree 
#using sample id
ID.tip.tree <- ggtree(tree) +
  geom_tiplab(size = 2.5, hjust = -0.1); ID.tip.tree
#ggsave(filename = "../244/output/ID_tip_tree.pdf", plot = ID.tip.tree, width = 20, height = 30, dpi = 1200)

#using host family
modify_tip(tree, "family")
family.tip.tree <- ggtree(tree) %<+% taxonomy.tip +
  geom_tiplab(size = 2, hjust = -0.1, aes(label = taxon)); family.tip.tree
#ggsave(filename = "../244/output/family_tip_tree.pdf", plot = family.tip.tree, width = 10, height = 15, dpi = 1200)

#using host taxa
modify_tip(tree, "mixed")
mixed.tip.tree <- ggtree(tree, ladderize = FALSE) %<+% taxonomy.tip +
  geom_treescale(x = 0, y = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +
  geom_tiplab(size = 4, hjust = -0.1, aes(label = taxon)) +
  geom_nodelab(aes(label = label), size = 2, hjust = -0.3); mixed.tip.tree
#ggsave(filename = "../244/output/mixed_tip_tree.pdf", plot = mixed.tip.tree, width = 10, height = 15, dpi = 1200)


################### by orthogroup #######################
orthogroup <- dataset.clean %>%
  group_by(orthogroups_in_OG, OGs_in_orthogroup, new_orthogroup) %>%
  mutate(gene=paste(toString(unique(Preferred_name)), toString(unique(Description)))) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  select(new_orthogroup, n, curated_COG_category, OGs_in_orthogroup, gene, mean_position_in_orthogroup, mean_GC_in_orthogroup, sd_GC_in_orthogroup, mean_length_in_orthogroup) %>%
  unique()

orthogroup.sample <- dataset.clean %>%
  select(new_orthogroup, sample, GC) %>%
  pivot_wider(names_from = sample, values_from = GC) %>%
  full_join(orthogroup, by = join_by(new_orthogroup)) %>%
  select(new_orthogroup, curated_COG_category, OGs_in_orthogroup, gene, n, everything()) %>%
  as.data.frame

################# omega #############################
M0.omega <- read.csv("../input/omega_summary.out", header = T)

M0.omega <- M0.omega %>%
  filter(omega < 2) 

M0.omega <- left_join(M0.omega, orthogroup.sample[c("new_orthogroup", "mean_position_in_orthogroup", "mean_GC_in_orthogroup", "sd_GC_in_orthogroup")], by = c("new_orthogroup"))


#combine KEGG information
source("KEGGtable.R")
ko.categories <- ko.table %>%
  select(!Preferred_name) %>%
  separate_rows(new_orthogroup, sep = ", ") %>%
  mutate(new_orthogroup = as.numeric(new_orthogroup)) %>%
  left_join(M0.omega, by = "new_orthogroup")

ko.categories <- ko.categories %>%
  filter(pathway %in% pathway.list) %>%
  mutate(category2 = factor(category2, levels = unique(category2))) %>%
  as.data.frame()
