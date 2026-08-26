syn <- read.table("../input/CAM19_1_C054_CR_S54_L001/output/tables/ms12syn.tsv", header = T)
nonsyn <- read.table("../input/CAM19_1_C054_CR_S54_L001/output/tables/ms12nonsyn.tsv", header = T)
all <- read.table("../input/CAM19_1_C054_CR_S54_L001/output/tables/ms12all.tsv", header = T)
combine <- full_join(syn, nonsyn, join_by(Mut)) 

combine <- combine %>%
  mutate(ObsNumAll = ObsNum.x + ObsNum.y) %>%
  mutate(ExpNumAll = ExpNum.x + ExpNum.y)

