source("omega/OmegaMain.R")

#GC
p1 <- ggplot(data = M0.omega, aes(x = mean_GC_in_orthogroup)) +
  theme_classic() +
  geom_density() + 
  xlab("GC-content") +
  ylab("Density") +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme(text = element_text(size = 14)); p1
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/GC1.png", plot = p1, width = 6, height = 4, dpi = 1200)
mean(M0.omega$mean_GC_in_orthogroup)
sd(M0.omega$mean_GC_in_orthogroup)
range(M0.omega$mean_GC_in_orthogroup)


p2 <- ggplot(data = M0.omega, aes(x = mean_position_in_orthogroup, y = mean_GC_in_orthogroup)) +
  theme_classic() +
  geom_point() + 
  xlab("Chromosomal position") +
  ylab("GC-content") +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme(text = element_text(size = 14)); p2
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/GC2.png", plot = p2, width = 6, height = 4, dpi = 1200)
cor.test(M0.omega$mean_GC_in_orthogroup, M0.omega$mean_position_in_orthogroup, method = "spearman")

p3 <- ggplot(data = M0.omega.COG, aes(y = mean_GC_in_orthogroup, x = curated_COG_category)) +
  theme_classic() +
  geom_boxplot() +
  xlab("COG category") +
  ylab("GC-content") +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega.COG)), hjust = 1.5, vjust = 2) +
  theme(text = element_text(size = 14)); p3
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/GC3.png", plot = p3, width = 6, height = 4, dpi = 1200)
kruskal.test(M0.omega$omega ~ M0.omega$curated_COG_category)

p4 <- ggplot(ko.categories, aes(x = category2, y = mean_GC_in_orthogroup)) +
  geom_boxplot(fill = "skyblue", color = "darkblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "blue") +  
  theme_classic() +
  labs(x = "KEGG categories", y = "GC-content") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  theme(text = element_text(size = 14)) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(ko.categories)), hjust = 1.5, vjust = 2); p4
#ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/GC4.png", plot = p4, width = 9, height = 6, dpi = 1200)


p5 <- ggplot(ko.categories, aes(x = factor(pathway, levels = unique(pathway)), y = mean_GC_in_orthogroup, fill = as.factor(category2))) +
  geom_boxplot(color = "darkblue", outlier.shape = NA) + 
  geom_jitter(width = 0.1, alpha = 0.6, color = "blue") + 
  theme_classic() +
  labs(x = "KEGG pathway", y = "GC-content", fill = "Category") +
  theme(axis.text.x = element_text(angle = 80, hjust = 1)) +
  theme(text = element_text(size = 14)) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(ko.categories)), hjust = 1.5, vjust = 2); p5
ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/GC5.png", plot = p5, width = 12, height = 8, dpi = 1200)
kruskal.test(ko.categories$mean_GC_in_orthogroup ~ ko.categories$pathway)

#test how much GC content is predicted by the gene and the sample
GC.content <- dataset.clean %>%
  select(query, GC, new_orthogroup, sample) %>%
  group_by(sample) %>%
  mutate(mean_GC_in_sample = mean(GC)) %>%
  ungroup() %>%
  group_by(new_orthogroup) %>%
  mutate(mean_GC_in_orthogroup = mean(GC)) %>%
  ungroup() %>%
  rename(orthogroup = new_orthogroup)


model.GC <- lmer(GC ~ (1 | sample) + (1 | orthogroup), data = GC.content)
summary(model.GC)
r2(model.GC)

tab_model(model.GC, 
          transform = NULL, 
          show.re.var = TRUE,
          linebreak = TRUE,
          dv.labels = c("Model: lmer(GC ~ (1 | sample) + (1 | orthogroup))"),
          pred.labels = c("(Intercept)"),
          file = "../244/output/geneByGene/geneGC.html")

png(filename = "~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/GC_random_effect.png",
  width = 1600,
  height = 1600,
  res = 300)
plot(GC.content$GC ~ fitted(model.GC),
     xlab = "Fitted GC-content",
     ylab = "Observed GC-content",
     pch = 16, cex = 0.4)
abline(0, 1, col = "red", lwd = 2)
legend("bottomright", legend = paste0("N = ", length(GC.content$GC)), bty = "n")
dev.off()








