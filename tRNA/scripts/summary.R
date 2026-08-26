source("main.R")

all_aa <- unique(tRNA_all$aa); all_aa
#20 amino acids (plus ile2 and fMet)

aa.anticodon <- tRNA_all %>%
  select(aa, anticodon) %>%
  unique()
#each amino acid has 1-4 anticodons


sample.tRNA <- tRNA_all %>%
  group_by(sample) %>%
  mutate(n = n()) %>%
  select(sample, n, family_combined) %>%
  unique()
#each sample has 31-34 tRNAs

all_tRNAs <- unique(tRNA_all$tRNA)
#34 unique tRNAs

per.tRNA.sample <- tRNA_all %>%
  group_by(tRNA, sample) %>%
  mutate(n = n())
#every gene (aa * anticodon combination) only has maximum one copy in each sample

per.tRNA <- tRNA_all %>%
  group_by(tRNA) %>%
  mutate(no.of.sample = n()) %>%
  select(tRNA, no.of.sample) %>%
  unique() %>%
  arrange(tRNA)
#most genes are conserved across all species, except: 
write.table(
  per.tRNA,
  "~/Desktop/Thesis/Typst/figures/chap3results/tRNAs.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE)

#find out the missing tRNAs
full_grid <- CJ(
  sample = unique(tRNA_all$sample),
  tRNA   = all_tRNAs)

missing_tRNAs <- full_grid[
  !tRNA_all[, .(sample, tRNA)],
  on = .(sample, tRNA)]

missing_tRNAs <- missing_tRNAs %>%
  left_join(taxonomy, by = join_by(sample == Sample))

  