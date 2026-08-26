library(tibble)

source("Main.R")
source("pathways/HeatmapFunction.R")

OGs <- c("3IV77@39782|Blattabacteriaceae", "3IVFI@39782|Blattabacteriaceae", "3IVG7@39782|Blattabacteriaceae", "3IVCE@39782|Blattabacteriaceae", "3IVGF@39782|Blattabacteriaceae", "3IVGS@39782|Blattabacteriaceae")
trpF <- "597"
trpB <- "598"


TrpPath <- orthogroup.sample %>%
  filter(OGs_in_orthogroup %in% OGs  | new_orthogroup == trpF | new_orthogroup == trpB) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column(var = "sample")

colnames(TrpPath) <- TrpPath[1, ]

TrpPath <- TrpPath %>%
#  filter(new_orthogroup %in% collapsed.tree$tip.label) %>%
#  filter(!complete.cases(.)) %>%
  mutate(family = taxonomy$family[match(new_orthogroup, taxonomy$sample)])

#gene gc heatmap
trp.operon <- c(426, 291, 418, 396, 597, 598, 115)
gc.heatmap(trp.operon)
