
########################### genes without annotation and containing only repeats ###########################################################
#CDS:CAM19_2_C100_CR_S4_L004_Anaplectidae_I03530 that is just repeats
condition <- grepl("CDS:CAM19_2_C100_CR_S4_L004", dataset$query) & is.na(dataset$eggNOG_OGs) & dataset$nucleotides > 500
dataset$notes[condition] <- "only repeats"

#CDS:CAM19_2_C171_CR_S70_L004_Anaplectidae_I03560 that is just repeats
condition <- grepl("CDS:CAM19_2_C171_CR_S70_L004", dataset$query) & is.na(dataset$eggNOG_OGs) & dataset$nucleotides > 500
dataset$notes[condition] <- "only repeats"

#CDS:AU035_S364_Pseudophyllodromiidae_I03190 repeats
condition <- grepl("CDS:AU035_S364", dataset$query) & is.na(dataset$eggNOG_OGs) & dataset$nucleotides > 500
dataset$notes[condition] <- "only repeats"

#CDS:AU069_S398_Pseudophyllodromiidae_I04650 repeats
condition <- grepl("CDS:AU069_S398", dataset$query) & is.na(dataset$eggNOG_OGs) & dataset$nucleotides > 500
dataset$notes[condition] <- "only repeats"

#CDS:ZU017_S88_L004_Blaberidae_I01290 repeats
condition <- grepl("CDS:ZU017_S88_L004", dataset$query) & is.na(dataset$eggNOG_OGs) & dataset$nucleotides > 500
dataset$notes[condition] <- "only repeats"

#CDS:AU006_S335_Pseudophyllodromiidae_I02200 repeats
condition <- grepl("CDS:AU006_S335", dataset$query) & is.na(dataset$eggNOG_OGs) & dataset$nucleotides > 600
dataset$notes[condition] <- "only repeats"

#CDS:AU033_S362_Pseudophyllodromiidae_I01040 repeats
condition <- grepl("CDS:AU033_S362", dataset$query) & is.na(dataset$eggNOG_OGs) & dataset$nucleotides > 800
dataset$notes[condition] <- "only repeats"

#CDS:CAM19_2_C005_CR_S7_L002_undetermined_I03620 repeats
condition <- grepl("CDS:CAM19_2_C005_CR_S7_L002", dataset$query) & is.na(dataset$eggNOG_OGs) & dataset$nucleotides > 800
dataset$notes[condition] <- "only repeats"

#CDS:CAM19_2_C020_CR_S21_L002_undetermined_I02020 repeats
condition <- grepl("CDS:CAM19_2_C020_CR_S21_L002", dataset$query) & is.na(dataset$eggNOG_OGs) & dataset$nucleotides > 600 & dataset$position > 200
dataset$notes[condition] <- "only repeats"

#CDS:BL940_S55_L003_Pseudophyllodromiidae_I02010 repeats
condition <- grepl("CDS:BL940_S55_L003", dataset$query) & is.na(dataset$eggNOG_OGs) & dataset$nucleotides > 1000
dataset$notes[condition] <- "only repeats"

#CDS:CAM19_2_C065_CR_S66_L002_Pseudophyllodromiidae_I02060 spanning a long gap
condition <- grepl("CDS:CAM19_2_C065_CR_S66_L002", dataset$query) & is.na(dataset$eggNOG_OGs) & dataset$nucleotides > 500
dataset$notes[condition] <- "only repeats"


########################### repeats within a gene or near the termini  ########################
#Gene rho/nusA have repeats within the gene in some samples
condition <- grepl("COG1158", dataset$eggNOG_OGs) & dataset$relative_length_in_orthogroup > 1.1
dataset$notes[condition] <- "contains repeats within the gene"

#GrpE some samples having repeats in the N-terminal disordered regions
condition <- grepl("COG0576", dataset$eggNOG_OGs)
dataset$notes[condition] <- "some samples have repeats near the N-terminus"

#mrcA some samples have repeats in the C-terminus 
condition <- grepl("COG5009", dataset$eggNOG_OGs)
dataset$notes[condition] <- "some samples have repeats near the C-terminus"

#CDS:CAM19_1_C019_CR_S19_L001_Anaplectidae_I03390 have a repeat in the C-terminus
condition <- grepl("CAM19_1_C019_CR_S19_L001", dataset$query) & grepl("COG0203", dataset$eggNOG_OGs)
dataset$notes[condition] <- "have repeats near the C-terminus"

#rpsS some samples have repeats in the C-terminus 
condition <- grepl("COG0185", dataset$eggNOG_OGs)
dataset$notes[condition] <- "some samples have repeats near the C-terminus"

#CDS:AU006_S335_Pseudophyllodromiidae_I01010 have a repeat in the C-terminus
condition <- grepl("AU006_S335", dataset$query) & grepl("COG0227", dataset$eggNOG_OGs)
dataset$notes[condition] <- "have repeats near the C-terminus"



#many more...



