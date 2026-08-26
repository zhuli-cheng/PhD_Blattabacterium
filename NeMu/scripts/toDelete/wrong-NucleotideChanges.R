source("Main.R")

nucleotide.changes <- data.frame(
  group = character(),
  A_net_gain = numeric(),
  T_net_gain = numeric(),
  G_net_gain = numeric(),
  C_net_gain = numeric(),
  stringsAsFactors = FALSE
)

for (group in selected.groups$group) {
  # Load mutation spectrums
  ms12syn <- read.table(paste0("../input/monophyleticGroups/NeMuPipeline/", group, "/output/tables/ms12syn.tsv"), header = TRUE)
  ms12all <- read.table(paste0("../input/monophyleticGroups/NeMuPipeline/", group, "/output/tables/ms12all.tsv"), header = TRUE)
  
  # Calculate the net change of the 4 nucleotides
  ms12all <- ms12all %>%
    separate(Mut, into = c("OriginalNucleotide", "DestinationNucleotide"), sep = ">")
  
  A.net.gain <- sum(ms12all$ObsNum[ms12all$DestinationNucleotide == "A"]) - sum(ms12all$ObsNum[ms12all$OriginalNucleotide == "A"])
  T.net.gain <- sum(ms12all$ObsNum[ms12all$DestinationNucleotide == "T"]) - sum(ms12all$ObsNum[ms12all$OriginalNucleotide == "T"])
  G.net.gain <- sum(ms12all$ObsNum[ms12all$DestinationNucleotide == "G"]) - sum(ms12all$ObsNum[ms12all$OriginalNucleotide == "G"])
  C.net.gain <- sum(ms12all$ObsNum[ms12all$DestinationNucleotide == "C"]) - sum(ms12all$ObsNum[ms12all$OriginalNucleotide == "C"])
  
  # Add a row to the results table
  nucleotide.changes <- rbind(nucleotide.changes, data.frame(
    group = group,
    A_net_gain = A.net.gain,
    T_net_gain = T.net.gain,
    G_net_gain = G.net.gain,
    C_net_gain = C.net.gain
  ))
}


