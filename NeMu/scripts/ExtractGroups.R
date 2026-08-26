source("LoadPackages.R")

tree <- read.tree("~/Desktop/BlattabacteriumRscripts/OrthologsAnalyses/244/input/244_RY.treefile")
two.families <- c("CAM19_1_C019_CR_S19_L001", "YS-3-1_S60_L003")
common.ancestor <- getMRCA(tree, two.families); common.ancestor
tree <- root(tree, node = common.ancestor, resolve.root = T)

extract_groups <- function(samples) {
  clade <- getMRCA(tree, samples)
  tips_from_clade <- extract.clade(tree, clade)
  tips <- tips_from_clade$tip.label
  return(tips)
}

#most groups are obtained by visually selecting monophyletic clusters from the Blattabacterium core protein tree of 244 samples;
#sometimes outgroups are distantly related when difficulty selecting outgroups within the clade
#some samples are from the mitochondrial phylogeny, when the Blattabacterium tree has too few samples, and the samples are verified to contain circular Blattabacterium genomes (among the 316 circular samples)


#Blattellidae1
samples <- c("CAM19_2_C091_CR_S89_L002", "CAM19_1_C073_CR_S72_L001")
group1 <- extract_groups(samples)
outgroup.group1 <- "CAM19_1_C073_CR_S72_L001"

#Blattellidae2 
samples <- c("OKE024_S37_L003", "CAM19_1_C010_CR_S10_L001")
group2 <- extract_groups(samples)
outgroup.group2 <- "CAM19_1_C010_CR_S10_L001"

#Blattellidae3
samples <- c("CAM19_1_C033_CR_S33_L001", "CAM19_1_C059_CR_S59_L001")
group3 <- extract_groups(samples)
outgroup.group3 <- "AU080_S409"
group3 <- c(group3, "AU080_S409")

#Blattellidae4
samples <- c("AU080_S409", "ZU021_S92_L004")
group4 <- extract_groups(samples)
outgroup.group4 <- "ZU021_S92_L004"

#Blaberidae1
samples <- c("ZU005_S76_L004", "ZU096_S8_L003")
group5 <- extract_groups(samples)
outgroup.group5 <- "ZU096_S8_L003"

#Blaberidae2
samples <- c("OKE026_S39_L003", "ZU091_S7_L003")
group6 <- extract_groups(samples)
outgroup.group6 <- "ZU096_S8_L003"
group6 <- c(group6, "ZU096_S8_L003")

#Blaberidae3
samples <- c("ZU013_S84_L004", "OKE010_S23_L003")
group7 <- extract_groups(samples)
outgroup.group7 <- "ZU091_S7_L003"
group7 <- c(group7, "ZU091_S7_L003")

#Blaberidae4
samples <- c("CAM19_2_C119_CR_S22_L004", "CAM19_1_C069_CR_S69_L001")
group8 <- extract_groups(samples)
outgroup.group8 <- "CAM19_2_C103_CR_S7_L004"
group8 <- c(group8, "CAM19_2_C103_CR_S7_L004")

#Blaberidae5
samples <- c("CAM19_2_C079_CR_S80_L002", "CAM19_2_C103_CR_S7_L004")
group9 <- extract_groups(samples)
outgroup.group9 <- "CAM19_2_C103_CR_S7_L004"

#Pseudophyllodromiidae1
samples <- c("AU006_S335", "AU057_S386")
group10 <- extract_groups(samples)
outgroup.group10 <- "AU008_S337" 
group10 <- c(group10, "AU008_S337")

#Pseudophyllodromiidae2
samples <- c("AU073_S402", "AU078_S407")
group11 <- extract_groups(samples)
outgroup.group11 <- "AU063_S392" 
group11 <- c(group11, "AU063_S392")

#Pseudophyllodromiidae3
group12 <- c("OKE002_S15_L003", "AU045_S374", "OKE023_S36_L003", "OKE028_S41_L003", "OKE022_S35_L003")
outgroup.group12 <- "ZU008_S79_L004" 
group12 <- c(group12, "ZU008_S79_L004")

#Pseudophyllodromiidae4
group12 <- c("OKE005_S18_L003", "CAM19_2_C065_CR_S66_L002", "CAM19_2_C090_CR_S88_L002", "CAM19_2_C018_CR_S19_L002", "YS-1fin_S58_L003", "OKE017_S30_L003", "CAM19_2_C020_CR_S21_L002")
outgroup.group12 <- "CAM19_1_C057_CR_S57_L001" 
group12 <- c(group12, "CAM19_1_C057_CR_S57_L001")

#Ectobiidae
samples <- c("AU074_S403", "BL937_S52_L003")
group13 <- extract_groups(samples)
outgroup.group13 <- "AU063_S392" 


#Blattidae
samples <- c("ZU016_S87_L004", "ZU019_S90_L004")
group14 <- extract_groups(samples)
outgroup.group14 <- "CAM19_2_C166_CR_S65_L004" 
group14 <- c(group14, "CAM19_2_C166_CR_S65_L004")


#Anaplectidae1
samples <- c("CAM19_2_C100_CR_S4_L004", "CAM19_1_C050_CR_S50_L001")
group15 <- extract_groups(samples)
outgroup.group15 <- "CAM19_1_C019_CR_S19_L001" 
group15 <- c(group15, "CAM19_1_C019_CR_S19_L001")

#Anaplectidae2
group16 <- c("AU082_S411", "AU076_S405", "AU087_S416", "BL928_S44_L003")
outgroup.group16 <- "CAM19_1_C019_CR_S19_L001" 
group16 <- c(group16, "CAM19_1_C019_CR_S19_L001")


#Corydiidae
group17 <- c("CAM19_1_C061_CR_S61_L001", "CAM19_1_C060_CR_S60_L001", "CAM19_1_C056_CR_S56_L001", "CAM19_2_C166_CR_S65_L004", 
             "CAM19_1_C041_CR_S41_L001", "YS-4_S61_L003", "YS-3-1_S60_L003", "YS-5_S62_L003")
outgroup.group17 <- "ZU016_S87_L004" 
group17 <- c(group17, "ZU016_S87_L004")

for (i in 1:17) {
  obj_name <- paste0("group", i)
  vec <- get(obj_name)
  writeLines(as.character(vec), paste0("../output/", obj_name, ".txt"))
}







