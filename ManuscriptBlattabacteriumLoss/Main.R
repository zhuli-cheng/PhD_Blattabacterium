source("LoadPackages.R")

source("../../TaxonTips/TaxonTips.R")

####################### read in input data #################################
sample.size <- 244

#read in full gene table
dataset.clean <- read.csv(paste0("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsCuration/", sample.size, "/output/dataset_clean.csv"), header = T)

#read in sample table
sample.orthogroup <- read.csv(paste0("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsCuration/", sample.size, "/output/sample_orthogroup.csv"), header = T)
sample.orthogroup[is.na(sample.orthogroup)] <- "0"
sample.orthogroup$family <- taxonomy$family[match(sample.orthogroup$sample, taxonomy$sample)]
sample.summary <- sample.orthogroup %>%
  select(sample, no_of_orthologs, mean_GC_of_orthologs)

#read in orthogroup table
orthogroup.sample <- read.csv(paste0("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsCuration/", sample.size, "/output/orthogroup_sample.csv"), header = T, check.names = FALSE)

#read in phylogeny 
tree <- read.tree("~/Desktop/BlattabacteriumRscripts/MitochondrialTree/input/MT106andRefSeqRY.treefile")

#reroot the tree (using Corydiidae and Anaplectidae)
two.families <- c("Anaplecta_omei_NC_070177_1", "Eucorydia_linglong_NC_082499_1")
common.ancestor <- getMRCA(tree, two.families); common.ancestor
tree <- root(tree, node = common.ancestor, resolve.root = T)
tree$node.label[tree$node.label == "Root"] <- ""


#keep only tips with complete Blattabacterium genomes
collapsed.tree <- drop.tip(tree, setdiff(tree$tip.label, sample.orthogroup$sample))

#process orthogroup.sample 
highly.similar <-  setdiff(tree$tip.label, collapsed.tree$tip.label)
orthogroup.sample <- orthogroup.sample[, !(names(orthogroup.sample) %in% highly.similar)]
orthogroup.sample$n <- rowSums(!is.na(orthogroup.sample[, grepl("_S", names(orthogroup.sample))]))

sample.orthogroup <- sample.orthogroup %>%
  filter(sample %in% collapsed.tree$tip.label)


#visualize the tree 
#using sample id
ID.tip.tree <- ggtree(collapsed.tree, ladderize = FALSE) +
  geom_tiplab(size = 2.5, hjust = -0.1); ID.tip.tree

#using host family
modify_tip(collapsed.tree, "family")
family.tip.tree <- ggtree(collapsed.tree) %<+% taxonomy.tip +
  geom_tiplab(size = 2, hjust = -0.1, aes(label = taxon)); family.tip.tree

#using host taxa
modify_tip(collapsed.tree, "mixed")
mixed.tip.tree <- ggtree(collapsed.tree, ladderize = FALSE) %<+% taxonomy.tip +
  geom_treescale(x = 0, y = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +
  geom_tiplab(size = 4, hjust = -0.1, aes(label = taxon)) +
  geom_nodelab(aes(label = label), size = 2, hjust = -0.3); mixed.tip.tree

