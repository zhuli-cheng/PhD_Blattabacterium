source("main.R")

ncRNA <- ncRNA_all %>%
  filter(!target_name %in% c("LSU_rRNA_bacteria", "SSU_rRNA_bacteria", "5S_rRNA", "tRNA")) %>%
  filter(E_value < 1.0e-5)

unique(ncRNA$target_name)

#find out their GA (Gathering Threshold) values from Rfam.cm
GA <- data.frame(
  target_name = c(
    "5_ureB_sRNA",
    "RNaseP_bact_a",
    "tmRNA",
    "Bacteria_small_SRP",
    "sau-50",
    "mir-593",
    "mir-574",
    "IRES_L-myc",
    "MIR1122",
    "Alfamo_CPB"
  ),
  GA = c(
    40,
    100,
    47,
    40,
    50,
    80,
    80,
    60,
    77,
    70
  ),
  stringsAsFactors = FALSE
)

#filter hits by GA
ncRNA <- ncRNA %>%
  left_join(GA, by = "target_name") %>%
           filter(score >= GA)

ncRNA <- ncRNA %>%
  group_by(sample) %>%
  mutate(n.RNA = n()) %>%
  ungroup() %>%
  group_by(target_name) %>%
  mutate(n.sample = n())

ncRNA.genes <- ncRNA %>%
  select(target_name, accession, n.sample) %>%
  unique() %>%
  filter(target_name != "sau-50")

write.table(
  ncRNA.genes,
  "~/Desktop/Thesis/Typst/figures/chap3results/ncRNAs.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE)




