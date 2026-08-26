add.ORF <- read.csv(paste0("../", sample.size, "/input/add_ORF.csv"), header = F)
colnames(add.ORF) <- c("new_orthogroup", "sample", "nucleotides", "GC")

add.ORF <- add.ORF %>%
  mutate(type = "ManuallyAdded") %>%
  mutate(location = "chromosome") %>%
  mutate(query = paste0("added:", sample, "_", new_orthogroup))
         

new_rows <- anti_join(add.ORF, dataset.clean, by = c("sample", "new_orthogroup"))
dataset.clean <- bind_rows(dataset.clean, new_rows) 
dataset.clean <- dataset.clean %>%
  group_by(new_orthogroup) %>%
  mutate(curated_COG_category = toString(unique(na.omit(curated_COG_category)))) %>%
  ungroup() 

dataset.clean <- summarize_orthogroups(dataset.clean, new_orthogroup)

#remove trpFB (fused in most samples) and add trpB and trpF
dataset.clean <- dataset.clean %>%
  filter(new_orthogroup != "437")

dataset.clean$curated_COG_category[dataset.clean$new_orthogroup == "597"] <- "E"
dataset.clean$Preferred_name[dataset.clean$new_orthogroup == "597"] <- "trpF"
dataset.clean$EC[dataset.clean$new_orthogroup == "597"] <- "4.2.1.20,5.3.1.24"
dataset.clean$KEGG_ko[dataset.clean$new_orthogroup == "597"] <- "ko:K01696,ko:K01817"
dataset.clean$KEGG_Pathway[dataset.clean$new_orthogroup == "597"] <- "ko00260,ko00400,ko01100,ko01110,ko01130,ko01230,map00260,map00400,map01100,map01110,map01130,map01230"
dataset.clean$KEGG_Module[dataset.clean$new_orthogroup == "597"] <- "ko00260,ko00400,ko01100,ko01110,ko01130,ko01230,map00260,map00400,map01100,map01110,map01130,map01230"
dataset.clean$mean_position_in_orthogroup[dataset.clean$new_orthogroup == "597"] <- "210.2"


dataset.clean$curated_COG_category[dataset.clean$new_orthogroup == "598"] <- "E"
dataset.clean$Preferred_name[dataset.clean$new_orthogroup == "598"] <- "trpB"
dataset.clean$EC[dataset.clean$new_orthogroup == "598"] <- "4.2.1.20,5.3.1.24"
dataset.clean$KEGG_ko[dataset.clean$new_orthogroup == "598"] <- "ko:K01696,ko:K01817"
dataset.clean$KEGG_Pathway[dataset.clean$new_orthogroup == "598"] <- "ko00260,ko00400,ko01100,ko01110,ko01130,ko01230,map00260,map00400,map01100,map01110,map01130,map01230"
dataset.clean$KEGG_Module[dataset.clean$new_orthogroup == "598"] <- "ko00260,ko00400,ko01100,ko01110,ko01130,ko01230,map00260,map00400,map01100,map01110,map01130,map01230"
dataset.clean$mean_position_in_orthogroup[dataset.clean$new_orthogroup == "598"] <- "210.1"
