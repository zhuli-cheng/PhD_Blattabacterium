############load ortholog data
setwd("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsAnalyses/scripts/")
source("Main.R")

############load taxonomy data
setwd("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/MitochondrialTree/scripts/")
source("../../TaxonTips/TaxonTips.R")

source("LoadPackages.R")

#############load tree and reroot
mito.tree <- read.tree("../input/MT106andRefSeqRY.treefile")
two.families <- c("Anaplecta_omei_NC_070177_1", "Eucorydia_linglong_NC_082499_1")
common.ancestor <- getMRCA(mito.tree, two.families); common.ancestor
mito.tree <- root(mito.tree, node = common.ancestor, resolve.root = T)
#mito.tree$node.label[mito.tree$node.label == "Root"] <- ""
ggtree(mito.tree, ladderize = FALSE) +
  geom_tiplab(size = 3, hjust = -0.1) 
  

#############load sequencing results on symbionts and gc content
absence <- read.csv("../input/blattabacterium_absence.csv")
colnames(absence) <- c("sample", "blattabacterium_absence", "blattabacterium", "wolbachia", "rickettsia", "absence_confirmed")
absence <- absence %>%
  mutate(wolbachia_presence = ifelse(wolbachia > 0, 1, 0)) %>%
  mutate(rickettsia_presence = ifelse(rickettsia > 0, 1, 0)) %>%
  mutate(rickettsiales_presence = ifelse(wolbachia > 0, 1, ifelse(rickettsia > 0, 1, 0)))

GC.content <- read.csv("../input/GCcontent.out", header = T)

GC.absence <- full_join(absence, GC.content, by = "sample")
rownames(GC.absence) <- GC.absence$sample
GC.absence <- GC.absence %>%
  select(!sample)
GC.absence[]<- lapply(GC.absence, as.numeric)
GC.absence[GC.absence == 0] <- NA

################modify tip labels
modify_tip(mito.tree, "mixed")
absence <- absence %>% 
  full_join(taxonomy.tip, join_by(sample == tip_label)) %>%
  mutate(taxon = ifelse(is.na(taxon) | taxon == "", sample, taxon))

absence.mixed <- absence %>%
  mutate(absence_confirmed = ifelse(is.na(absence_confirmed), 0, absence_confirmed)) %>%
  mutate(absence_confirmed = as.factor(absence_confirmed))

absence.mixed$absence_confirmed[absence.mixed$sample == "CAM19_2_C164_CR_S63_L004"] <- "1"
absence.mixed$absence_confirmed[absence.mixed$sample == "CAM19_1_C037_CR_S37_L001"] <- "1"
absence.mixed$absence_confirmed[absence.mixed$sample == "CAM19_2_C007_CR_S9_L002"] <- "1"

#######################plot tree
ID.tip.tree <- ggtree(mito.tree, ladderize = FALSE) %<+% absence.mixed +
  geom_tiplab(size = 3, hjust = -0.1, aes(label = taxon, color = absence_confirmed)) +
  geom_treescale(x = 0, y = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +
  geom_text2(aes(subset = !isTip, label = label),
             hjust = -0.2,
             vjust = 0.5,
             size = 2.5) +
  scale_color_manual(values = c("0" = "black", "1"  = "blue")) + 
  theme(legend.position = "none"); ID.tip.tree
#ggsave(filename = "../output/BlattabacteriumLossTreeSmall.pdf", plot = ID.tip.tree, width = 10, height = 16, dpi = 1200)
#ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/BlattabacteriumLossTreeSmall.pdf", plot = ID.tip.tree, width = 10, height = 14, dpi = 1200)

########################3other symbionts
tips_to_remove <- grep("NC_", mito.tree$tip.label, value = TRUE)
mito.tree <- drop.tip(mito.tree, tips_to_remove)

symbionts <- absence.mixed %>%
  select(sample, absence_confirmed, wolbachia_presence, rickettsia_presence)
rownames(symbionts) <- symbionts$sample
symbionts$sample <- NULL

symbionts <- symbionts[mito.tree$tip.label, ,]
symbionts[is.na(symbionts)] <- 0
symbionts$absence_confirmed <- as.numeric(as.character(symbionts$absence_confirmed))
symbionts <- symbionts %>%
  mutate(blattabacterium_presence = ifelse(absence_confirmed == 0, 2, 0)) %>%
  select(blattabacterium_presence, wolbachia_presence, rickettsia_presence)

svg("../output/symbionts.svg", width = 8, height = 6)
phylo.heatmap(mito.tree, symbionts, colors = 0:7, fsize = 0.01, standardize = FALSE, legend = FALSE, labels = FALSE, split=c(0.8, 0.2))
dev.off()

#writeLines(mito.tree$tip.label, "../output/tip_labels.txt")

reads.446 <- read.csv("../input/Assembly_Bacteria_reads.csv")
reads.filtered <- reads.446[reads.446$sample %in% rownames(symbionts), ]
reads.filtered <- reads.filtered %>%
  select(sample, Blattabacterium_proportion, B_depth, Wolbachia_proportion, W_depth, Rickettsia_proportion, R_depth)

write.table(
  reads.filtered,
  "~/Desktop/Thesis/Typst/figures/supplementary materials/Assembly_Bacteria_reads_106.csv",
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = TRUE)



