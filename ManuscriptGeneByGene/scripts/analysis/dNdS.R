source("omega/OmegaMain.R")

#remove the gene (amiD) with omega > 1
M0.omega <- M0.omega %>%
  filter(omega < 1)

#omega
p.omega1 <- ggplot(data = M0.omega, aes(x = omega)) +
  theme_classic() +
  geom_density() + 
  xlab(expression(italic(d)[N] / italic(d)[S])) +
  ylab("Density") +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme_AnimalBehaviour(); p.omega1
mean(M0.omega$omega)
sd(M0.omega$omega)
range(M0.omega$omega)
sum(M0.omega$omega < 0.1)
sum(M0.omega$omega < 0.2)


p.omega2 <- ggplot(data = M0.omega, aes(x = mean_position_in_orthogroup, y = omega)) +
  theme_classic() +
  geom_point() + 
  xlab("Chromosomal position") +
  ylab(expression(italic(d)[N] / italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme_AnimalBehaviour(); p.omega2
cor.test(M0.omega$omega, M0.omega$mean_position_in_orthogroup, method = "spearman")


p4 <- ggplot(ko.categories, aes(x = category2, y = omega)) +
  geom_boxplot(fill = "skyblue", color = "darkblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "blue") +  
  theme_classic() +
  labs(x = "KEGG categories", y = expression(italic(d)[N] / italic(d)[S] * " of the orthogroup")) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  theme(text = element_text(size = 14)) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(ko.categories)), hjust = 1.5, vjust = 2); p4
kruskal.test(ko.categories$omega ~ ko.categories$category2)


p.omega5 <- ggplot(ko.categories, aes(x = factor(pathway, levels = unique(pathway)), y = omega)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, alpha = 0.6) + 
  theme_classic() +
  labs(x = "KEGG pathway", y = expression(italic(d)[N] / italic(d)[S]), fill = "Category") +
#  theme(axis.text.x = element_text(angle = 80, hjust = 1)) +
  theme_AnimalBehaviour() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = "none") +
  theme(text = element_text(size = 14)) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(ko.categories)), hjust = 1.5, vjust = 2); p.omega5


kruskal.test(ko.categories$omega ~ ko.categories$pathway)


