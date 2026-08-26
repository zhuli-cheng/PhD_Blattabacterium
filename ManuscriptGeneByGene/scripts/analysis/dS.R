source("omega/OmegaMain.R")

#remove the gene (amiD) with omega > 1
M0.omega <- M0.omega %>%
  filter(omega < 1)

ko.categories <- ko.categories %>%
  filter(omega < 1)

#omega
p.dS1 <- ggplot(data = M0.omega, aes(x = dS_standardized)) +
  theme_classic() +
  geom_density() + 
  xlab(expression("Standardized " * italic(d)[S])) +
  ylab("Density") +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme_AnimalBehaviour(); p.dS1
mean(M0.omega$dS_standardized)
sd(M0.omega$dS_standardized)
range(M0.omega$dS_standardized)

p.dS2 <- ggplot(data = M0.omega, aes(x = mean_position_in_orthogroup, y = dS_standardized)) +
  theme_classic() +
  geom_point() + 
  xlab("Chromosomal position") +
  ylab(expression("Standardized " * italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  geom_vline(xintercept = 323, linetype = "dashed", size = 1) +
  annotate("text", x = 323, y = 130, label = "rRNA operon", vjust = 1.5, angle = 90, size = 5) +
  theme_AnimalBehaviour(); p.dS2



p4 <- ggplot(ko.categories, aes(x = category2, y = dS_standardized)) +
  geom_boxplot(fill = "skyblue", color = "darkblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "blue") +  
  theme_classic() +
  labs(x = "KEGG categories", y = expression("Orthogroup " * italic(d)[S] * " (relative to core protein tree length)")) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  theme(text = element_text(size = 14)) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(ko.categories)), hjust = 1.5, vjust = 2); p4
kruskal.test(ko.categories$dS_standardized ~ ko.categories$category2)


p.dS5 <- ggplot(ko.categories, aes(x = factor(pathway, levels = unique(pathway)), y = dS_standardized, fill = as.factor(category2))) +
  geom_boxplot(color = "darkblue", outlier.shape = NA) + 
  geom_jitter(width = 0.1, alpha = 0.6, color = "blue") + 
  theme_classic() +
  labs(x = "KEGG category", y = expression("Standardized " * italic(d)[S])) +
  theme(axis.text.x = element_text(angle = 80, hjust = 1)) +
  theme(text = element_text(size = 14)) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(ko.categories)), hjust = 1.5, vjust = 2); p.dS5
kruskal.test(ko.categories$dS_standardized ~ ko.categories$pathway)
