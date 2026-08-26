#this manual curation includes the CDS and PPgenes from the 245 new genome (circular and very few gaps)


########################### modifications on the annotations/OGs etc ###########################################################

dataset$curated_OG[dataset$OG_Blattabacteriaceae == "3IVMI@39782|Blattabacteriaceae"] <- "3IVHQ@39782|Blattabacteriaceae"

dataset$curated_OG[dataset$OG_Blattabacteriaceae == "3IVKC@39782|Blattabacteriaceae"] <- "3IVKW@39782|Blattabacteriaceae"

dataset$curated_OG[dataset$OG_Flavobacteriia == "1HX0Q@117743|Flavobacteriia"] <- "1IJ6G@117743|Flavobacteriia"
dataset$curated_OG[dataset$OG_Flavobacteriia == "1IJ6G@117743|Flavobacteriia"] <- "1IJ6G@117743|Flavobacteriia"

dataset$curated_OG[dataset$OG_root == "COG0198@1|root"] <- "3IVKY@39782|Blattabacteriaceae"

dataset$curated_OG[dataset$OG_Flavobacteriia == "1HYI1@117743|Flavobacteriia"] <- "1HYI1@117743|Flavobacteriia"

#dataset$curated_COG_category[dataset$OG_Blattabacteriaceae == "3IVKD@39782|Blattabacteriaceae"] <- "S"

#dataset$curated_COG_category[dataset$OG_Blattabacteriaceae == "3IVHY@39782|Blattabacteriaceae"] <- "S"

dataset$curated_COG_category[dataset$OG_Blattabacteriaceae == "3IVKI@39782|Blattabacteriaceae"] <- "E"

#dataset$curated_COG_category[dataset$OG_Blattabacteriaceae == "3IVJ6@39782|Blattabacteriaceae"] <- "S"

#dataset$curated_COG_category[dataset$OG_Blattabacteriaceae == "3IVJS@39782|Blattabacteriaceae"] <- "S"

dataset$curated_COG_category[dataset$OG_Blattabacteriaceae == "3IVKV@39782|Blattabacteriaceae"] <- "C"

dataset$curated_COG_category[dataset$OG_Blattabacteriaceae == "3IVJ6@39782|Blattabacteriaceae"] <- "J"





########################### delete samples ###########################
#CAM19_1_C074_CR_S73_L001_undetermined_I03380, problem with the ends of the original contig. And many redundant genomes. already deleted before curation
#dataset <- dataset %>%
#  filter(!(grepl("CAM19_1_C074_CR_S73_L001", query)))

#sample_size <- sample_size - 1


################################ take notes on possible duplication/error###########################
#CDS:YS-1fin_S58_L003_undetermined_I03520 possible duplication or error (next to a gap)
duplicated <- dataset %>%
  filter((grepl("CDS:YS-1fin_S58_L003", query) & grepl("COG0094", eggNOG_OGs))) %>%
  arrange(position) %>%
  slice(n())
dataset <- anti_join(dataset, duplicated, by = join_by(query))
duplicated$notes <- "possible duplication or error and next to a gap"
dataset <- merge(dataset, duplicated, all = TRUE)

#CDS:CAM19_1_C014_CR_S14_L001_Blaberidae_I03940 possible duplication or error
duplicated <- dataset %>%
  filter((grepl("CDS:CAM19_1_C014_CR_S14_L001", query) & grepl("COG0089", eggNOG_OGs))) %>%
  arrange(position) %>%
  slice(n())
dataset <- anti_join(dataset, duplicated, by = join_by(query))
duplicated$notes <- "possible duplication or error and next to a gap"
dataset <- merge(dataset, duplicated, all = TRUE)

#PPgene:ZU091_S7_L003_Blaberidae_I01280 possible duplication of a short region or error
condition <- grepl("ZU091_S7_L003", dataset$query) & grepl("COG0158", dataset$eggNOG_OGs) 
dataset$notes[condition] <- "duplication of a short region or error within the PPgene"

#CDS:AU045_S374_Pseudophyllodromiidae_I04020 and CDS:AU045_S374_Pseudophyllodromiidae_I04030 possible duplication (same aa, slightly different nt). no gaps nearby
condition <- grepl("AU045_S374", dataset$query) & grepl ("2ACTF", dataset$eggNOG_OGs)
duplicated <- dataset[condition, ] %>%
  arrange(position) %>%
  slice(n())
dataset <- anti_join(dataset, duplicated, by = join_by(query))
duplicated$notes <- "possible duplication"
dataset <- merge(dataset, duplicated, all = TRUE)

#CDS:AU062_S391_Ectobiidae_I03860 error due to gap
condition <- grepl("CDS:AU062_S391", dataset$query) & is.na(dataset$eggNOG_OGs)
error <- dataset[condition, ] %>%
  subset(position > 300) %>%
  subset(position < 400)
query <- error$query
condition <- grepl(query, dataset$query)
dataset$notes[condition] <- "false duplication: error due to a gap"
dataset <- dataset[!condition,]

#CDS:OKE016_S29_L003_Blaberidae_I02220 error due to being the end of the assembly
condition <- grepl("OKE016_S29_L003", dataset$query) & grepl ("2ACTD", dataset$eggNOG_OGs)
error <- dataset[condition, ] %>%
  arrange(position) %>%
  slice(n())
query <- error$query
condition <- grepl(query, dataset$query)
dataset$notes[condition] <- "false duplication: error due to being end of the assembly"


############## delete the false CDS within rRNA ############
condition <- grepl("_I0", dataset$query) & dataset$GC > 50
dataset$notes[condition] <- "rRNA"


