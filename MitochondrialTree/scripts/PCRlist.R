source("Main.R")

PCR.list <- read.csv("../input/PCR.list", header = F)
colnames(PCR.list) <- "sample"

PCR.list.host <- PCR.list %>%
  left_join(taxonomy, by = "sample") %>%
  select(mixed)

write.table(
  PCR.list.host,
  "../output/PCR.list.host.csv",
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE)
