library(tibble)

source("Main.R")
source("Heatmaps.R")

OGs <- c("3IVGT@39782|Blattabacteriaceae", "3IV9I@39782|Blattabacteriaceae", "3IVDW@39782|Blattabacteriaceae", "3IV8U@39782|Blattabacteriaceae", "3IV42@39782|Blattabacteriaceae", "3IVCZ@39782|Blattabacteriaceae", "3IV5M@39782|Blattabacteriaceae", "3IV59@39782|Blattabacteriaceae")

HisPath <- orthogroup.sample %>%
  filter(OGs_in_orthogroup %in% OGs) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample")

colnames(HisPath) <- HisPath[1, ]

HisPath <- HisPath %>%
  #  filter(new_orthogroup %in% collapsed.tree$tip.label) %>%
  #  filter(!complete.cases(.)) %>%
  mutate(family = taxonomy$family[match(new_orthogroup, taxonomy$sample)])

#gene gc heatmap
his.operon <- c(310, 170, 341, 193, 61, 438, 8, 49)
gc.heatmap(his.operon)




