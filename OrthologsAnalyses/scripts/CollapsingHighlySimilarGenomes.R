#define function to collapse highly similar strains
collapse_clades <- function(samples) {
  clade <- getMRCA(collapsed.tree, samples)
  tips_from_clade <- extract.clade(collapsed.tree, clade)
  drop.tip(collapsed.tree, tips_from_clade$tip.label[2:length(tips_from_clade$tip.label)])
}

#Blattellidae1
samples <- c("CAM19_2_C056_CR_S57_L002", "CAM19_2_C156_CR_S55_L004")
collapsed.tree <- collapse_clades(samples)

#Blattellidae2
samples <- c("CAM19_2_C002_CR_S4_L002", "CAM19_2_C114_CR_S17_L004")
collapsed.tree <- collapse_clades(samples)

#Blattellidae3
samples <- c("OKE006_S19_L003", "OKE009_S22_L003")
collapsed.tree <- collapse_clades(samples)

#Blattellidae4
samples <- c("OKE020_S33_L003", "ZU072_S3_L003")
collapsed.tree <- collapse_clades(samples)

#Blattellidae5
samples <- c("ZU109_S10_L003", "OKE016_S29_L003")
collapsed.tree <- collapse_clades(samples)

#Blattellidae6
samples <- c("ZU018_S89_L004", "ZU004_S75_L004")
collapsed.tree <- collapse_clades(samples)

#Blattellidae7
samples <- c("CAM19_2_C170_CR_S69_L004", "ZU076_S4_L003")
collapsed.tree <- collapse_clades(samples)

#Blattellidae8
samples <- c("ZU059_S1_L003", "CAM19_2_C144_CR_S45_L004")
collapsed.tree <- collapse_clades(samples)

#Blattellidae9
samples <- c("CAM19_1_C058_CR_S58_L001", "CAM19_2_C169_CR_S68_L004")
collapsed.tree <- collapse_clades(samples)

#Blattellidae10
samples <- c("CAM19_2_C056_CR_S57_L002", "CAM19_2_C112_CR_S15_L004")
collapsed.tree <- collapse_clades(samples)

#Blattellidae11
samples <- c("ZU009_S80_L004", "ZU013_S84_L004")
collapsed.tree <- collapse_clades(samples)

#Blattellidae12
samples <- c("ZU109_S10_L003", "OKE020_S33_L003")
collapsed.tree <- collapse_clades(samples)

#Blattellidae13
samples <- c("CAM19_2_C170_CR_S69_L004", "CAM19_1_C069_CR_S69_L001")
collapsed.tree <- collapse_clades(samples)

#Blattellidae14
samples <- c("CAM19_1_C058_CR_S58_L001", "ZU059_S1_L003")
collapsed.tree <- collapse_clades(samples)





#Blaberidae1
samples <- c("CAM19_2_C121_CR_S24_L004", "CAM19_1_C089_CR_S86_L001")
collapsed.tree <- collapse_clades(samples)

#Blaberidae2
samples <- c("CAM19_1_C066_CR_S66_L001", "CAM19_1_C076_CR_S75_L001")
collapsed.tree <- collapse_clades(samples)

#Blaberidae3
samples <- c("CAM19_2_C122_CR_S25_L004", "CAM19_2_C024_CR_S25_L002")
collapsed.tree <- collapse_clades(samples)

#Blaberidae4
samples <- c("CAM19_2_C096_CR_S94_L002", "CAM19_1_C001_CR_S1_L001")
collapsed.tree <- collapse_clades(samples)

#Blaberidae5
samples <- c("CAM19_1_C068_CR_S68_L001", "CAM19_1_C072_CR_S71_L001")
collapsed.tree <- collapse_clades(samples)

#Blaberidae6
samples <- c("CAM19_2_C159_CR_S58_L004", "CAM19_2_C153_CR_S52_L004")
collapsed.tree <- collapse_clades(samples)

#Blaberidae7
samples <- c("CAM19_2_C064_CR_S65_L002", "CAM19_2_C075_CR_S76_L002")
collapsed.tree <- collapse_clades(samples)

#Blaberidae8
samples <- c("CAM19_2_C068_CR_S69_L002", "CAM19_2_C012_CR_S14_L002")
collapsed.tree <- collapse_clades(samples)

#Blaberidae9
samples <- c("CAM19_2_C029_CR_S30_L002", "CAM19_2_C046_CR_S47_L002")
collapsed.tree <- collapse_clades(samples)

#Blaberidae10
samples <- c("CAM19_1_C042_CR_S42_L001", "CAM19_2_C082_CR_S82_L002")
collapsed.tree <- collapse_clades(samples)

#Blaberidae11
samples <- c("CAM19_2_C051_CR_S52_L002", "CAM19_1_C043_CR_S43_L001")
collapsed.tree <- collapse_clades(samples)

#Blaberidae12
samples <- c("CAM19_1_C016_CR_S16_L001", "CAM19_1_C044_CR_S44_L001")
collapsed.tree <- collapse_clades(samples)

#Blaberidae13
samples <- c("AU049_S378", "AU054_S383")
collapsed.tree <- collapse_clades(samples)


#Blaberidae14
samples <- c("CAM19_1_C066_CR_S66_L001", "CAM19_1_C008_CR_S8_L001")
collapsed.tree <- collapse_clades(samples)

#Blaberidae15
samples <- c("CAM19_2_C096_CR_S94_L002", "CAM19_1_C086_CR_S84_L001")
collapsed.tree <- collapse_clades(samples)

#Blaberidae16
samples <- c("YS-2-1_S59_L003", "ZU001_S72_L004")
collapsed.tree <- collapse_clades(samples)

#Blaberidae17
samples <- c("CAM19_1_C016_CR_S16_L001", "CAM19_2_C030_CR_S31_L002")
collapsed.tree <- collapse_clades(samples)

#Blaberidae18
samples <- c("AU055_S384", "AU080_S409")
collapsed.tree <- collapse_clades(samples)



#Pseudophyllodromiidae1
samples <- c("OKE023_S36_L003", "OKE028_S41_L003")
collapsed.tree <- collapse_clades(samples)


#Pseudophyllodromiidae2
samples <- c("AU057_S386", "AU011_S340")
collapsed.tree <- collapse_clades(samples)

#Pseudophyllodromiidae3
samples <- c("OKE028_S41_L003", "OKE022_S35_L003")
collapsed.tree <- collapse_clades(samples)

#Ectobiidae1
samples <- c("AU027_S356", "AU018_S347")
collapsed.tree <- collapse_clades(samples)

#Ectobiidae2
samples <- c("AU002_S331", "AU044_S373")
collapsed.tree <- collapse_clades(samples)

#Anaplectidae1
samples <- c("CAM19_2_C100_CR_S4_L004", "CAM19_2_C171_CR_S70_L004")
collapsed.tree <- collapse_clades(samples)

#Blattidae1
samples <- c("CAM19_2_C133_CR_S35_L004", "CAM19_2_C125_CR_S28_L004")
collapsed.tree <- collapse_clades(samples)

#Blattidae2
samples <- c("ZU020_S91_L004", "ZU022_S93_L004")
collapsed.tree <- collapse_clades(samples)






ggtree(collapsed.tree) +
  geom_tiplab(size = 2, hjust = -0.1)
#compare AU080_S409 and AU055_S384
#selected.rows <- subset(ortholog.table, X %in% c("AU080_S409", "AU055_S384"))
#selected.rows <- subset(ortholog.table, X %in% c("AU096_S425", "AU002_S331", "AU062_S391"))
#selected.rows <- subset(ortholog.table, X %in% c("OKE022_S35_L003", "OKE023_S36_L003", "OKE028_S41_L003"))
#selected.rows <- subset(ortholog.table, X %in% c("CAM19_2_C171_CR_S70_L004", "CAM19_2_C100_CR_S4_L004"))
#selected.rows <- subset(ortholog.table, X %in% c("ZU112_S11_L003", "ZU020_S91_L004", "ZU022_S93_L004"))
#selected.rows <- as.data.frame(lapply(selected.rows, as.numeric))
#difference <- selected.rows[2, ] - selected.rows[1, ]
#difference <- difference[, -c(1, 2)]
#difference <- difference[, -ncol(difference)]
#different <- difference[, difference[1, ] != 0]

