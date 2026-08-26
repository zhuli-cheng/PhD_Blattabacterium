source("omega/OmegaMain.R")

#remove the gene (amiD) with omega > 1
M0.omega <- M0.omega %>%
  filter(omega < 1)
M0.omega.COG <- M0.omega.COG %>%
  filter(omega < 1)
ko.categories <- ko.categories %>%
  filter(omega < 1)

#omega
p1 <- ggplot(data = M0.omega, aes(x = omega)) +
  theme_classic() +
  geom_density() + 
  xlab(expression(italic(d)[N] / italic(d)[S])) +
  ylab("Density") +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme(text = element_text(size = 14)); p1
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/omega1.png", plot = p1, width = 6, height = 4, dpi = 1200)
mean(M0.omega$omega)
sd(M0.omega$omega)
range(M0.omega$omega)
sum(M0.omega$omega < 0.1)
sum(M0.omega$omega > 0.2)


p2 <- ggplot(data = M0.omega, aes(x = mean_position_in_orthogroup, y = omega)) +
  theme_classic() +
  geom_point() + 
  xlab("Chromosomal position") +
  ylab(expression(italic(d)[N] / italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme(text = element_text(size = 14)); p2
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/omega2.png", plot = p2, width = 6, height = 4, dpi = 1200)
cor.test(M0.omega$omega, M0.omega$mean_position_in_orthogroup, method = "spearman")


p3 <- ggplot(data = M0.omega.COG, aes(y = omega, x = curated_COG_category)) +
  theme_classic() +
  geom_boxplot() +
  xlab("COG category") +
  ylab(expression(italic(d)[N] / italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega.COG)), hjust = 1.5, vjust = 2) +
  theme(text = element_text(size = 14)); p3
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/omega3.png", plot = p3, width = 6, height = 4, dpi = 1200)
kruskal.test(M0.omega$omega ~ M0.omega$curated_COG_category)


p4 <- ggplot(ko.categories, aes(x = category2, y = omega)) +
  geom_boxplot(fill = "skyblue", color = "darkblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "blue") +  
  theme_classic() +
  labs(x = "KEGG categories", y = expression(italic(d)[N] / italic(d)[S] * " of the orthogroup")) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  theme(text = element_text(size = 14)) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(ko.categories)), hjust = 1.5, vjust = 2); p4
#ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/omega4.png", plot = p4, width = 9, height = 6, dpi = 1200)
kruskal.test(ko.categories$omega ~ ko.categories$category2)


p5 <- ggplot(ko.categories, aes(x = factor(pathway, levels = unique(pathway)), y = omega, fill = as.factor(category2))) +
  geom_boxplot(color = "darkblue", outlier.shape = NA) + 
  geom_jitter(width = 0.1, alpha = 0.6, color = "blue") + 
  theme_classic() +
  labs(x = "KEGG pathway", y = expression(italic(d)[N] / italic(d)[S]), fill = "Category") +
  theme(axis.text.x = element_text(angle = 80, hjust = 1)) +
  theme(text = element_text(size = 14)) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(ko.categories)), hjust = 1.5, vjust = 2); p5
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/omega5.png", plot = p5, width = 12, height = 8, dpi = 1200)
kruskal.test(ko.categories$omega ~ ko.categories$pathway)
