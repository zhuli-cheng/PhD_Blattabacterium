library(dplyr)
library(stringr)
library(ggplot2)
library(lme4)
library(tidyr)
library(tibble)

#choose different sample collections: 152 for 98% similarity threshold (old + new samples); 244 for new samples with minimum gaps; 383 for all circular genomes (old and new)
sample.size <- 244

########################### step 1. read data ###########################################################
#read information on host family
family <- read.csv("../../sample_family.csv")
colnames(family) <- c("sample", "host_family", "host_family_inferred", "species")
family <- family %>% 
  mutate(combined = ifelse(host_family_inferred != "", paste0(host_family_inferred, "Inferred"), host_family_inferred))
family$combined[family$combined == ""] <- family$host_family[family$combined == ""]
family$host_family <- as.character(family$combined)
family <- family %>%
  select(sample, host_family)


#read stats on gene length(nucleotides) and GC content
representatives <- readLines(paste0("../", sample.size, "/input/list"))

gene.stats <- read.table(paste0("../", sample.size, "/input/gene_stats.tsv"), header = F) %>%
  setNames(c("query", "nucleotides", "GC"))

Proteinortho <- read.table(paste0("../", sample.size,"/input/proteinortho.tsv"), header = F) %>%
  setNames(c("query", "orthogroup"))

#read annotations from eggNOG
eggNOG <- read.csv(paste0("../", sample.size, "/input/383_diamond.emapper.annotations"),
                 sep = "\t",
                 comment.char = "#",
                 header = F,
                 na.strings = "-")

colnames(eggNOG) <- c("query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
                    "max_annot_lvl", "COG_category", "Description", "Preferred_name",
                    "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction",
                    "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
                    "BiGG_Reaction", "PFAMs")


#combine proteinortho, eggNOG, gene stats and host family
dataset <- full_join(eggNOG, gene.stats, by = "query") %>%
  full_join(Proteinortho, by = "query") %>%
  mutate(sample = str_extract(query, "(?<=:).*?(?=_[^_]*$)")) %>%
  left_join(family, by = "sample")

#filter representative genomes
dataset <- dataset[dataset$sample %in% representatives, ]


########################### step 2. manual curation ###########################################################

dataset$notes <- NA

source("GeneFusion.R")
  
dataset$OG_root <- str_extract(dataset$eggNOG_OGs, "\\w+@1\\|root")
dataset$OG_Blattabacteriaceae <- str_extract(dataset$eggNOG_OGs, "\\w+@39782\\|Blattabacteriaceae")
dataset$OG_Flavobacteriia <- str_extract(dataset$eggNOG_OGs, "\\w+@117743\\|Flavobacteriia")

##add curated columns
dataset$curated_COG_category <- dataset$COG_category
dataset$curated_OG <- dataset$OG_Blattabacteriaceae
dataset$curated_orthogroup <- dataset$orthogroup

dataset <- dataset %>%
  mutate(position = as.integer(str_extract(dataset$query, "(?<=_I)[0-9]+"))/10)

source("VariousCuration.R")
source("IncompleteCDS.R")

########################### step 3. perform some dataset manipulations to obtain a one-to-one match between orthogroups and annotations ###########################################################
dataset <- dataset %>%
  mutate(type = str_extract(dataset$query, "CDS|PPgene")) %>%
  mutate(location = str_extract(dataset$query, "(?<=_)(?!.*_)[^_0-9]+(?=[0-9])")) %>%
  mutate(location = ifelse(location == "I", "chromosome", ifelse(location == "II", "plasmid", location))) 


dataset <- dataset %>%
  group_by(curated_orthogroup) %>%
  mutate(position_min = min(position)) %>%
  mutate(position_max = max(position)) %>%
  mutate(position_mean = mean(position)) %>%
  ungroup() 


#define a function that will calculate extra statistics 
summarize_orthogroups <- function(dataframe, choose_orthogroup) {
    dataframe %>%
    group_by({{ choose_orthogroup }}) %>%
    mutate(genes_in_orthogroup = n()) %>%
    mutate(sample_size_in_orthogroup = n_distinct(sample)) %>%
    mutate(mean_GC_in_orthogroup = mean(na.omit(GC))) %>%
    mutate(sd_GC_in_orthogroup = sd(na.omit(GC))) %>%
    mutate(mean_length_in_orthogroup = mean(na.omit(nucleotides))) %>%
    mutate(relative_length_in_orthogroup = nucleotides/mean_length_in_orthogroup) %>%
    mutate(distance_length_in_orthogroup = nucleotides - mean_length_in_orthogroup) %>%
    mutate(max_length_in_orthogroup = max(na.omit(relative_length_in_orthogroup))) %>%
    mutate(min_length_in_orthogroup = min(na.omit(relative_length_in_orthogroup))) %>%
    mutate(OGs_in_orthogroup = ifelse(length(unique(curated_OG[!is.na(curated_OG)])) == 0, NA, 
                                      toString(unique(curated_OG[!is.na(curated_OG)])))) %>% 
    mutate(OGs_in_orthogroup = ifelse(!is.na(OGs_in_orthogroup), OGs_in_orthogroup, ifelse(!is.na(unique(OG_root[!is.na(OG_root)])), toString(unique(OG_root[!is.na(OG_root)])), NA))) %>%
    mutate(locations_in_orthogroup = toString(unique(location))) %>%
    ungroup() %>%
    mutate(OG_inferred_from_orthogroup = ifelse(!is.na(curated_OG), curated_OG, OGs_in_orthogroup)) %>%
    group_by(OG_inferred_from_orthogroup) %>%
    mutate(orthogroups_in_OG = ifelse(is.na(OG_inferred_from_orthogroup), NA, toString(unique({{ choose_orthogroup }})))) %>%
    ungroup() %>%
    mutate(core_proteome = ifelse(sample_size_in_orthogroup == sample.size, TRUE, FALSE)) %>%
    group_by({{ choose_orthogroup }}) %>%
    mutate(curated_COG_category = toString(unique(na.omit(curated_COG_category)))) %>%
    ungroup() 
  }

#use the function on the dataset
dataset <- summarize_orthogroups(dataset, curated_orthogroup)

source("Repeats.R")

View(dataset[!is.na(dataset$notes),])


#recalculate statistics using the merged orthogroups(instead of the original orthogroups)
#only orthogroups not found on plasmid; only orthogroups with annotations
#delete pseudogenes
dataset.clean <- dataset %>%
  subset(!grepl("not the complete ortholog: only domains", notes))

dataset.clean <- dataset.clean %>%
  subset(!grepl("rRNA", notes))

dataset.clean <- dataset.clean %>%
  subset(locations_in_orthogroup == "chromosome")

dataset.clean <- summarize_orthogroups(dataset.clean[complete.cases(dataset.clean$OG_inferred_from_orthogroup), ], orthogroups_in_OG)


#with these statistics, check one-to-one relationship between orthogroups(merged) and OGs
orthogroup.OG <- dataset.clean %>%
  count(OGs_in_orthogroup, orthogroups_in_OG)

orthogroup.OG$OGs_in_orthogroup[duplicated(orthogroup.OG$OGs_in_orthogroup)]
orthogroup.OG$orthogroups_in_OG[duplicated(orthogroup.OG$orthogroups_in_OG)]

########################### step 4. look at the abnormal genes/orthogroups ###########################################################
##### genes and orthogroups with no annotations
no.annotation <- dataset %>%
  subset(is.na(OG_inferred_from_orthogroup))  %>%
  subset(locations_in_orthogroup == "chromosome") %>%
  subset(is.na(notes))

View(no.annotation %>%
  count(orthogroup, sample_size_in_orthogroup, mean_length_in_orthogroup, location))

#write.csv(no.annotation,"no.annotation.csv")

############ singletons
singletons <- dataset.clean %>%
  subset(genes_in_orthogroup == 1)

####### paralogs
putative.paralogs <- dataset.clean %>%
  group_by(orthogroups_in_OG) %>%
  filter(duplicated(sample) | duplicated(sample, fromLast = TRUE)) %>%
  mutate(paralogs = toString(query)) %>%
  ungroup() 

#delete sequence with the lowest scores
to.delete <- putative.paralogs %>%
  group_by(sample) %>%
  slice_min(score)

dataset.clean <- dataset.clean %>%
  subset(!query %in% to.delete$query)

#reorder the orthogroups
dataset.clean <- dataset.clean %>%
  group_by(OGs_in_orthogroup) %>%
  mutate(new_orthogroup = cur_group_id()) %>%
  ungroup()

length(unique(dataset.clean$orthogroups_in_OG))
length(unique(dataset.clean$OGs_in_orthogroup))
range(dataset.clean$new_orthogroup)


############################## save datasets #######################

#1. generate tables
source("GenerateTablesForAnalyses.R")
#2. manually curate in Bash (with blast, alignment and orf predictions)
#3. add ORFs
source("AddORF.R")
#4. generate tables again
source("GenerateTablesForAnalyses.R")


























######################## check gene length  ########################

###########look at gene length distribution
plot(density(dataset.clean$relative_length_in_orthogroup[!is.na(dataset.clean$relative_length_in_orthogroup)]))

#check gene length distribution in each orthogroup
for (group in 99:1) {
  data.subset <- filter(dataset.clean, new_orthogroup == group)
  plot <- ggplot(data.subset, aes(x = nucleotides/3)) +
    geom_boxplot() +
    xlim(0, NA) +
    ggtitle(paste("Density plot for", group))
  print(plot)
}

long.genes <- dataset.clean %>%
  subset(distance_length_in_orthogroup > 180) %>%
  subset(is.na(notes))
  
short.genes <- dataset.clean %>%
  subset(distance_length_in_orthogroup < -180) %>%
  subset(is.na(notes))

#write.csv(long.genes, "long.genes.csv")
#look at sequences shorter than expected

condition <- dataset.clean$distance_length_in_orthogroup < -180 | (dataset.clean$distance_length_in_orthogroup < -120 & dataset.clean$relative_length_in_orthogroup < 0.7)
dataset.clean$pseudogenized <- ifelse(condition, TRUE, FALSE)

condition <- dataset.clean$pseudogenized == FALSE & dataset.clean$distance_length_in_orthogroup < -120
dataset.clean$pseudogenized[condition] <- "Ambiguous"


short.genes <- dataset.clean %>%
  subset(pseudogenized == FALSE) %>%
  subset(distance_length_in_orthogroup > -120) %>%
  subset(distance_length_in_orthogroup < -90) %>%
#  subset(relative_length_in_orthogroup < 0.7) %>%
#  subset(min_length_in_orthogroup > 0.6) %>%
  subset(is.na(notes))

unique(short.genes$orthogroups_in_OG)

ggplot(short.genes, aes(x = curated_COG_category)) +
  geom_bar()

plot(density(na.omit(dataset.clean$distance_length_in_orthogroup)))


