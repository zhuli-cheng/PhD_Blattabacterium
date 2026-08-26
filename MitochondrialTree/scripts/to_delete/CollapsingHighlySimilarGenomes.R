#only use when using mito.tree <- read.tree("../input/iqtree.RY/MT448andRefSeqRY.treefile")

#define function to collapse highly similar strains
collapse_clades <- function(samples) {
  clade <- getMRCA(mito.tree, samples)
  tips_from_clade <- extract.clade(mito.tree, clade)
  drop.tip(mito.tree, tips_from_clade$tip.label[2:length(tips_from_clade$tip.label)])
}

samples <- c("CAM19_1_C005_CR_S5_L001","CAM19_1_C079_CR_S77_L001")
mito.tree <- collapse_clades(samples)






