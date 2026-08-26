source("Main.R")

#GC
p.GC1 <- ggplot(data = M0.omega, aes(x = mean_GC_in_orthogroup)) +
  theme_classic() +
  geom_density() + 
  xlab("GC content") +
  ylab("Density") +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme_AnimalBehaviour(); p.GC1
mean(M0.omega$mean_GC_in_orthogroup)
sd(M0.omega$mean_GC_in_orthogroup)
range(M0.omega$mean_GC_in_orthogroup)


p.GC2 <- ggplot(data = M0.omega, aes(x = mean_position_in_orthogroup, y = mean_GC_in_orthogroup)) +
  theme_classic() +
  geom_point() + 
  xlab("Chromosomal position") +
  ylab("GC content") +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme_AnimalBehaviour(); p.GC2
cor.test(M0.omega$mean_GC_in_orthogroup, M0.omega$mean_position_in_orthogroup, method = "spearman")


p4 <- ggplot(ko.categories, aes(x = category2, y = mean_GC_in_orthogroup)) +
  geom_boxplot(fill = "skyblue", color = "darkblue", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "blue") +  
  theme_classic() +
  labs(x = "KEGG categories", y = "GC-content") +
  theme(axis.text.x = element_text(angle = 70, hjust = 1)) +
  theme(text = element_text(size = 14)) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(ko.categories)), hjust = 1.5, vjust = 2); p4


p.GC5 <- ggplot(ko.categories, aes(x = factor(pathway, levels = unique(pathway)), y = mean_GC_in_orthogroup)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(width = 0.1, alpha = 0.6) + 
  labs(x = "KEGG pathway", y = "GC content", fill = "KEGG Category") +
#  theme(axis.text.x = element_text(angle = 80, hjust = 1)) +
  theme(text = element_text(size = 14)) +
  theme_AnimalBehaviour() +
  theme(axis.text.x = element_blank()) +
  theme(legend.position = "none") +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(ko.categories)), hjust = 1.5, vjust = 2); p.GC5
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
          file = "../output/geneGC.html")

png(filename = "../output/GC_random_effect.png",
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








