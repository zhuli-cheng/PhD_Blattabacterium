source("omega/OmegaMain.R")

#remove the gene (amiD) with omega > 1
M0.omega <- M0.omega %>%
  filter(omega < 1)
M0.omega.COG <- M0.omega.COG %>%
  filter(omega < 1)
ko.categories <- ko.categories %>%
  filter(omega < 1)

#omega
p1 <- ggplot(data = M0.omega, aes(x = dS_standardized)) +
  theme_classic() +
  geom_density() + 
  xlab(expression("Standardized " * italic(d)[S])) +
  ylab("Density") +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme(text = element_text(size = 14)); p1
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/dS1.png", plot = p1, width = 6, height = 4, dpi = 1200)
mean(M0.omega$dS_standardized)
sd(M0.omega$dS_standardized)
range(M0.omega$dS_standardized)

p2 <- ggplot(data = M0.omega, aes(x = mean_position_in_orthogroup, y = dS_standardized)) +
  theme_classic() +
  geom_point() + 
  xlab("Chromosomal position") +
  ylab(expression("Standardized " * italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  geom_vline(xintercept = 323, linetype = "dashed", size = 1) +
  annotate("text", x = 323, y = 40, label = "rDNA operon", vjust = 1.5, angle = 90, size = 5) +
  theme(text = element_text(size = 14)); p2
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/dS2.png", plot = p2, width = 6, height = 4, dpi = 1200)


p3 <- ggplot(data = M0.omega.COG, aes(y = dS_standardized, x = curated_COG_category)) +
  theme_classic() +
  geom_boxplot() +
  xlab("COG category") +
  ylab(expression("Standardized " * italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega.COG)), hjust = 1.5, vjust = 2) +
  theme(text = element_text(size = 14)); p3
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/dS3.png", plot = p3, width = 6, height = 4, dpi = 1200)
kruskal.test(M0.omega$dS_standardized ~ M0.omega$curated_COG_category)


p4 <- ggplot(ko.categories, aes(x = category2, y = dS_standardized)) +
  geom_boxplot(fill = "skyblue", color = "darkblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "blue") +  
  theme_classic() +
  labs(x = "KEGG categories", y = expression("Orthogroup " * italic(d)[S] * " (relative to core protein tree length)")) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  theme(text = element_text(size = 14)) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(ko.categories)), hjust = 1.5, vjust = 2); p4
#ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/dS4.png", plot = p4, width = 9, height = 6, dpi = 1200)
kruskal.test(ko.categories$dS_standardized ~ ko.categories$category2)


p5 <- ggplot(ko.categories, aes(x = factor(pathway, levels = unique(pathway)), y = dS_standardized, fill = as.factor(category2))) +
  geom_boxplot(color = "darkblue", outlier.shape = NA) + 
  geom_jitter(width = 0.1, alpha = 0.6, color = "blue") + 
  theme_classic() +
  labs(x = "KEGG category", y = expression("Standardized " * italic(d)[S])) +
  theme(axis.text.x = element_text(angle = 80, hjust = 1)) +
  theme(text = element_text(size = 14)) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(ko.categories)), hjust = 1.5, vjust = 2); p5
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/dS5.png", plot = p5, width = 12, height = 8, dpi = 1200)
kruskal.test(ko.categories$dS_standardized ~ ko.categories$pathway)

#dS rate heterogeneity
M0.dS <- M0.omega %>%
  filter(mean_position_in_orthogroup > 250 & mean_position_in_orthogroup < 400) %>%
  mutate(group = ifelse(mean_position_in_orthogroup > 323, "upstream", "downstream"))
#the breaking point is where the rRNA operon is

p6 <- ggplot(data = M0.dS, aes(x = mean_position_in_orthogroup, y = dS_standardized, color = group)) +
  geom_point() +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  geom_vline(xintercept = 323, linetype = "dashed", size = 1) +
  annotate("text", x = 323, y = 40, label = "rDNA operon", vjust = 1.5, angle = 90, size = 5) +
  labs(x = "Chromosomal position", y = expression("Standardized " * italic(d)[S]), color = "Region"); p6
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/dS6.png", plot = p6, width = 12, height = 8, dpi = 1200)

wilcox.test(dS_standardized ~ group, data = M0.dS,  alternative = "less")
