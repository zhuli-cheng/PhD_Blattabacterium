########################### gene fusion ################################################################################
#tyrA, pheB CDS:CAM19_2_C065_CR_S66_L002_Pseudophyllodromiidae_I02620
condition <- grepl("CAM19_2_C065_CR_S66_L002", dataset$query) & grepl("COG0287", dataset$eggNOG_OGs)

fused.gene <- dataset[condition,]
dataset <- dataset %>%
  filter(!condition)

fused.gene[, !names(fused.gene) %in% c("query", "GC", "nucleotides", "sample")] <- NA
fused.gene[, "eggNOG_OGs"] <- "COG0287@1|root,COG0287@2|Bacteria,4NEKF@976|Bacteroidetes,1HWRY@117743|Flavobacteriia,3IV8T@39782|Blattabacteriaceae"
fused.gene$notes <- "gene fusion"
orthogroup <- min(na.omit(dataset[grepl("COG0287", dataset$eggNOG_OGs), "orthogroup"]))
fused.gene[, "orthogroup"] <- orthogroup
dataset <- rbind(dataset, fused.gene)

fused.gene[, !names(fused.gene) %in% c("query", "GC", "nucleotides", "sample")] <- NA
fused.gene[, "eggNOG_OGs"] <- "COG1605@1|root,COG2876@1|root,COG1605@2|Bacteria,COG2876@2|Bacteria,4NDU4@976|Bacteroidetes,1HX6I@117743|Flavobacteriia,3IVB8@39782|Blattabacteriaceae"
fused.gene$notes <- "gene fusion"
orthogroup <- min(na.omit(dataset[grepl("COG1605", dataset$eggNOG_OGs), "orthogroup"]))
fused.gene[, "orthogroup"] <- orthogroup
dataset <- rbind(dataset, fused.gene)

#pssA, dnaQ CDS:AU045_S374_Pseudophyllodromiidae_I00960
condition <- grepl("AU045_S374", dataset$query) & grepl("COG0847", dataset$eggNOG_OGs)

fused.gene <- dataset[condition,]
dataset <- dataset %>%
  filter(!condition)

fused.gene[, !names(fused.gene) %in% c("query", "GC", "nucleotides", "sample")] <- NA
fused.gene[, "eggNOG_OGs"] <- "COG0847@1|root,COG0847@2|Bacteria,4NE82@976|Bacteroidetes,1HXW3@117743|Flavobacteriia,3IVCP@39782|Blattabacteriaceae"
fused.gene$notes <- "gene fusion"
orthogroup <- min(na.omit(dataset[grepl("COG0847", dataset$eggNOG_OGs), "orthogroup"]))
fused.gene[, "orthogroup"] <- orthogroup
dataset <- rbind(dataset, fused.gene)

fused.gene[, !names(fused.gene) %in% c("query", "GC", "nucleotides", "sample")] <- NA
fused.gene[, "eggNOG_OGs"] <- "COG1183@1|root,COG1183@2|Bacteria,4NNUZ@976|Bacteroidetes,1HZ4M@117743|Flavobacteriia,3IV44@39782|Blattabacteriaceae"
fused.gene$notes <- "gene fusion"
orthogroup <- min(na.omit(dataset[grepl("COG1183", dataset$eggNOG_OGs), "orthogroup"]))
fused.gene[, "orthogroup"] <- orthogroup
dataset <- rbind(dataset, fused.gene)

# false gene fusion caused by Ns; rpsA, nrdA CDS:BL940_S55_L003_Pseudophyllodromiidae_I02760
condition <- grepl("BL940_S55_L003", dataset$query) & grepl("COG0209", dataset$eggNOG_OGs)

fused.gene <- dataset[condition,]
dataset <- dataset %>%
  filter(!condition)

fused.gene[, !names(fused.gene) %in% c("query", "GC", "nucleotides", "sample")] <- NA
fused.gene[, "eggNOG_OGs"] <- "COG0209@1|root,COG0209@2|Bacteria,4NEHQ@976|Bacteroidetes,1HXHB@117743|Flavobacteriia,3IVAK@39782|Blattabacteriaceae"
fused.gene$notes <- "false gene fusion"
orthogroup <- min(na.omit(dataset[grepl("COG0209", dataset$eggNOG_OGs), "orthogroup"]))
fused.gene[, "orthogroup"] <- orthogroup
dataset <- rbind(dataset, fused.gene)

fused.gene[, !names(fused.gene) %in% c("query", "GC", "nucleotides", "sample")] <- NA
fused.gene[, "eggNOG_OGs"] <- "COG0539@1|root,COG0539@2|Bacteria,4NDW9@976|Bacteroidetes,1HXPG@117743|Flavobacteriia,3IV8X@39782|Blattabacteriaceae"
fused.gene$notes <- "false gene fusion"
orthogroup <- min(na.omit(dataset[grepl("COG0539", dataset$eggNOG_OGs), "orthogroup"]))
fused.gene[, "orthogroup"] <- orthogroup
dataset <- rbind(dataset, fused.gene)
