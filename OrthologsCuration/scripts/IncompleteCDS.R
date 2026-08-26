########################### take notes on some domains of potentially PPgenes to Domain #########################
#TamB gene, having only a few domains found in certain Blaberidae samples
condition <- grepl("COG2911", dataset$eggNOG_OGs) & dataset$nucleotides < 3000
dataset$notes[condition] <- "not the complete ortholog: only domains"

#ilvK gene, some of the queries are domains
condition <- grepl("COG0683", dataset$eggNOG_OGs) & dataset$nucleotides < 360
queries <- dataset$sample[condition]
condition <- dataset$sample %in% queries & grepl("COG0683", dataset$eggNOG_OGs)
dataset$notes[condition] <- "not the complete ortholog: only domains"

#LptC CDS:CAM19_1_C010_CR_S10_L001_undetermined_I05320 some of the queries are domains
condition <- grepl("COG3117", dataset$eggNOG_OGs) & grepl("CAM19_1_C010_CR_S10_L001", dataset$query)
dataset$notes[condition] <- "PPgene"

#rnpA only complete in a few samples
condition <- grepl("COG0594", dataset$eggNOG_OGs) | grepl("3IVMG", dataset$eggNOG_OGs)
orthogroups <- unique(dataset$orthogroup[condition])
condition <- dataset$orthogroup %in% orthogroups & dataset$nucleotides < 300
dataset$notes[condition] <- "not the complete ortholog: only domains"

#CDS:BL940_S55_L003_Pseudophyllodromiidae_I04780 and CDS:BL940_S55_L003_Pseudophyllodromiidae_I04790 partial domains, too short
condition <- grepl("BL940_S55_L003", dataset$query) & grepl("COG0535", dataset$eggNOG_OGs)
dataset$notes[condition] <- "not the complete ortholog: only domains"


#CDS:AU063_S392_Pseudophyllodromiidae_I04920 lost most domains to repeats
condition <- grepl("AU063_S392", dataset$query) & grepl("3IV4T", dataset$eggNOG_OGs)
dataset$notes[condition] <- "not the complete ortholog: only domains; lost most domains to repeats"

#gene spr/MepS lost the signal peptide in some CDS
condition1 <- grepl("COG0791", dataset$eggNOG_OGs) & dataset$relative_length_in_orthogroup > 1.3
condition2 <- grepl("COG0791", dataset$eggNOG_OGs) & dataset$relative_length_in_orthogroup <= 1.3
dataset$notes[condition1] <- "signal and chain peptides"
dataset$notes[condition2] <- "only chain peptide without signal peptide"

#LptC-related losing signal peptide in some CDS
condition1 <- grepl("COG3117", dataset$eggNOG_OGs) & dataset$relative_length_in_orthogroup > 1.08
condition2 <- grepl("COG3117", dataset$eggNOG_OGs) & dataset$relative_length_in_orthogroup <= 1.08 & is.na(dataset$notes)
dataset$notes[condition1] <- "signal and chain peptides"
dataset$notes[condition2] <- "only chain peptide without signal peptide"



