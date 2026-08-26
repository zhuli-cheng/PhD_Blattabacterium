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
tree <- read.tree("../244/input/244_RY.treefile")

#reroot the tree (using Corydiidae and Anaplectidae)
two.families <- c("CAM19_1_C019_CR_S19_L001", "YS-3-1_S60_L003")
common.ancestor <- getMRCA(tree, two.families); common.ancestor
tree <- root(tree, node = common.ancestor, resolve.root = T)
tree$node.label[tree$node.label == "Root"] <- ""

#tree <- drop.tip(tree, "OKE011_S24_L003") #OKE011 has controversial position on the tree

#collapse highly similar genomes; run when needed
collapsed.tree <- tree
source("CollapsingHighlySimilarGenomes.R")

#tree <- drop.tip(collapsed.tree, tip = "OKE011_S24_L003")
#write.tree(collapsed.tree, file = "../244/output/244_RY_pruned.tree")

#process orthogroup.sample 
highly.similar <-  setdiff(tree$tip.label, collapsed.tree$tip.label)
orthogroup.sample <- orthogroup.sample[, !(names(orthogroup.sample) %in% highly.similar)]
orthogroup.sample$n <- rowSums(!is.na(orthogroup.sample[, grepl("_S", names(orthogroup.sample))]))

sample.orthogroup <- sample.orthogroup %>%
  filter(sample %in% collapsed.tree$tip.label)

#write.tree(collapsed.tree, file=paste0("../", sample.size, "/output/collapsed_tree.nwk"))

#visualize the tree 
#using sample id
ID.tip.tree <- ggtree(collapsed.tree) +
  geom_tiplab(size = 2.5, hjust = -0.1); ID.tip.tree
#ggsave(filename = "../244/output/ID_tip_tree.pdf", plot = ID.tip.tree, width = 20, height = 30, dpi = 1200)

#using host family
modify_tip(collapsed.tree, "family")
family.tip.tree <- ggtree(collapsed.tree) %<+% taxonomy.tip +
  geom_tiplab(size = 2, hjust = -0.1, aes(label = taxon)); family.tip.tree
#ggsave(filename = "../244/output/family_tip_tree.pdf", plot = family.tip.tree, width = 10, height = 15, dpi = 1200)

#using host taxa
modify_tip(collapsed.tree, "mixed")
mixed.tip.tree <- ggtree(collapsed.tree, ladderize = FALSE) %<+% taxonomy.tip +
  geom_treescale(x = 0, y = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +
  geom_tiplab(size = 4, hjust = -0.1, aes(label = taxon)) +
  geom_nodelab(aes(label = label), size = 2, hjust = -0.3); mixed.tip.tree
#ggsave(filename = "../244/output/mixed_tip_tree.pdf", plot = mixed.tip.tree, width = 10, height = 15, dpi = 1200)

ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/Blattabacterium_PCGs_mixed_tip_tree.pdf", plot = mixed.tip.tree, width = 10, height = 14, dpi = 1200)

