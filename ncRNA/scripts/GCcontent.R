source("main.R")

rna23S <- ncRNA_all %>%
  filter(target_name == "LSU_rRNA_bacteria")
#gc 0.41-0.45

rna16S <- ncRNA_all %>%
  filter(target_name == "SSU_rRNA_bacteria")
#gc 0.45-0.48

rna5S <- ncRNA_all %>%
  filter(target_name == "5S_rRNA")
#gc 0.45-0.48

tRNA <- ncRNA_all %>%
  filter(target_name == "tRNA")
#gc 0.34-0.65
