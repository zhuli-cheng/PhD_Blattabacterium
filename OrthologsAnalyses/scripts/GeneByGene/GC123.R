source("omega/OmegaMain.R")
M0.omega <- M0.omega %>% 
  filter(omega < 1)

#load GC data
ffn <- list.files(path = "../244/input/ffn", pattern = "\\.ffn$", full.names = TRUE)

GC.all <- lapply(ffn, function(f) {
  fa <- read.fasta(f)
  data.frame(orthogroup = sub("\\.ffn$", "", basename(f)), contig = names(fa), GC  = sapply(fa, GC), GC1 = sapply(fa, GC1), GC2 = sapply(fa, GC2), GC3 = sapply(fa, GC3))
})

GC.all <- do.call(rbind, GC.all)

rownames(GC.all) <- NULL

GC.orthogroup <- GC.all %>%
  group_by(orthogroup) %>%
  mutate(mean_GC_of_orthogroup = mean(GC)) %>%
  mutate(mean_GC1_of_orthogroup = mean(GC1)) %>%
  mutate(mean_GC2_of_orthogroup = mean(GC2)) %>%
  mutate(mean_GC3_of_orthogroup = mean(GC3)) %>%
  ungroup() %>%
  select(orthogroup, mean_GC_of_orthogroup, mean_GC1_of_orthogroup, mean_GC2_of_orthogroup, mean_GC3_of_orthogroup) %>%
  unique()
  
ggplot() +
  theme_classic() +
  geom_density(data = GC.orthogroup, aes(x = mean_GC1_of_orthogroup, colour = "GC1"), linewidth = 1) +
  geom_density(data = GC.orthogroup, aes(x = mean_GC2_of_orthogroup, colour = "GC2"), linewidth = 1) +
  geom_density(data = GC.orthogroup, aes(x = mean_GC3_of_orthogroup, colour = "GC3"), linewidth = 1) +
  xlab("GC-content") +
  ylab("Density") +
  labs(colour = "Codon position") +
  theme(text = element_text(size = 14))

#combine with omega
GC.orthogroup$orthogroup <- as.integer(GC.orthogroup$orthogroup)
GC.omega <- full_join(M0.omega, GC.orthogroup, join_by(new_orthogroup == orthogroup))

plot(GC.omega$mean_GC1_of_orthogroup + GC.omega$mean_GC2_of_orthogroup, GC.omega$mean_GC3_of_orthogroup)

plot(GC.omega$omega, GC.omega$mean_GC1_of_orthogroup)
plot(GC.omega$omega, GC.omega$mean_GC2_of_orthogroup)
plot(GC.omega$omega, GC.omega$mean_GC3_of_orthogroup)

#write.csv(results, "gc_by_codon_position.csv", row.names = FALSE)
#ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/GC1.png", plot = p1, width = 6, height = 4, dpi = 1200)
