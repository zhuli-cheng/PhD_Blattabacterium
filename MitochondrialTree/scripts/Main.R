

#load ortholog data
setwd("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsAnalyses/scripts/")
source("Main.R")

#load taxonomy data
setwd("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/MitochondrialTree/scripts/")
source("../../TaxonTips/TaxonTips.R")
source("LoadPackages.R")

#load mitochondrial tree
#mito.tree <- read.tree("../input/MT450.treefile")
#mito.tree <- read.tree("../input/MT448andRefSeq.treefile")
#mito.tree <- read.tree("../input/iqtree.RY/MT448andRefSeqRY.contree")
mito.tree <- read.tree("../input/MT448andRefSeqRY.treefile")


two.families <- c("CAM19_2_C060_CR_S61_L002", "CAM19_1_C061_CR_S61_L001")
common.ancestor <- getMRCA(mito.tree, two.families); common.ancestor
mito.tree <- root(mito.tree, node = common.ancestor, resolve.root = T)
mito.tree$node.label[mito.tree$node.label == "Root"] <- ""

mito.tree <- drop.tip(mito.tree, tip = c("CAM19_1_C046_CR_S46_L001"))
mito.tree <- drop.tip(mito.tree, tip = c("OKE022_S35_L003"))
#remove 2_C133 which has only one gene
mito.tree <- drop.tip(mito.tree, "CAM19_2_C133_CR_S35_L004")



#load sequencing results on symbionts
absence <- read.csv("../input/blattabacterium_absence.csv")
colnames(absence) <- c("sample", "blattabacterium_absence", "blattabacterium", "wolbachia", "rickettsia", "absence_confirmed")
absence <- absence %>%
  mutate(wolbachia_presence = ifelse(wolbachia > 0, 1, 0)) %>%
  mutate(rickettsia_presence = ifelse(rickettsia > 0, 1, 0)) %>%
  mutate(rickettsiales_presence = ifelse(wolbachia > 0, 1, ifelse(rickettsia > 0, 1, 0)))

#load GC content of Blattabacterium
GC.content <- read.csv("../input/GCcontent.out", header = T)

#combine all
GC.absence <- full_join(absence, GC.content, by = "sample")
rownames(GC.absence) <- GC.absence$sample
GC.absence <- GC.absence %>%
  select(!sample)
GC.absence[]<- lapply(GC.absence, as.numeric)
GC.absence[GC.absence == 0] <- NA




