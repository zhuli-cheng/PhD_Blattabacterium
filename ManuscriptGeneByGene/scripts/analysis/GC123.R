source("Main.R")
M0.omega <- M0.omega %>% 
  filter(omega < 1)

#load GC data
ffn <- list.files(path = "../input/ffn", pattern = "\\.ffn$", full.names = TRUE)

GC.all <- lapply(ffn, function(f) {
  fa <- read.fasta(f)
  data.frame(orthogroup = sub("\\.ffn$", "", sub("^orthogroup_", "", basename(f))), contig = names(fa), GC  = sapply(fa, GC), GC1 = sapply(fa, GC1), GC2 = sapply(fa, GC2), GC3 = sapply(fa, GC3))
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
GC.omega <- left_join(M0.omega, GC.orthogroup, join_by(new_orthogroup == orthogroup))

plot(GC.omega$mean_GC1_of_orthogroup + GC.omega$mean_GC2_of_orthogroup, GC.omega$mean_GC3_of_orthogroup)

p.correlation4 <- ggplot(data = GC.omega, aes(y = mean_GC1_of_orthogroup, x = omega)) +
  geom_point() +
  ylab(expression(atop("GC content", "1"^"st" * " codon position"))) +
  xlab(expression(italic(d)[N]/italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(GC.omega)), hjust = 1.5, vjust = 2) +
  theme_AnimalBehaviour(); p.correlation4

p.correlation5 <- ggplot(data = GC.omega, aes(y = mean_GC2_of_orthogroup, x = omega)) +
  geom_point() +
  ylab(expression(atop("GC content", "2"^"nd" * " codon position"))) +
  xlab(expression(italic(d)[N]/italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(GC.omega)), hjust = 1.5, vjust = 2) +
  theme_AnimalBehaviour(); p.correlation5

p.correlation6 <- ggplot(data = GC.omega, aes(y = mean_GC3_of_orthogroup, x = omega)) +
  geom_point() +
  ylab(expression(atop("GC content", "3"^"rd" * " codon position"))) +
  xlab(expression(italic(d)[N]/italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(GC.omega)), hjust = 1.5, vjust = 2) +
  theme_AnimalBehaviour(); p.correlation6





