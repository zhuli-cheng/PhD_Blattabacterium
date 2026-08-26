source("omega/omegaMain.R")

#remove the gene (amiD) with omega > 1
M0.omega <- M0.omega %>%
  filter(omega < 1)
M0.omega.COG <- M0.omega.COG %>%
  filter(omega < 1)
ko.categories <- ko.categories %>%
  filter(omega < 1)


############## correlations ################
p1 <- ggplot(data = M0.omega, aes(y = mean_GC_in_orthogroup, x = dN_standardized)) +
  geom_point() +
  ylab("GC-content") +
  xlab(expression("Standardized " * italic(d)[N])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme_set(theme_classic(base_size = 14)); p1

p2 <- ggplot(data = M0.omega, aes(y = mean_GC_in_orthogroup, x = dS_standardized)) +
  geom_point() +
  ylab("GC-content") +
  xlab(expression("Standardized " * italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme_set(theme_classic(base_size = 14)); p2

p3 <- ggplot(data = M0.omega, aes(y = mean_GC_in_orthogroup, x = omega)) +
  geom_point() +
  ylab("GC-content") +
  xlab(expression(italic(d)[N]/italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme_set(theme_classic(base_size = 14)); p3
cor.test(x = M0.omega$mean_GC_in_orthogroup, y = M0.omega$omega, method = "pearson")

p <- p1 / p2 / p3
p
ggsave("~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/GC_correlations.png", plot = p, width = 8, height = 12)

p4 <- ggplot(data = M0.omega, aes(x = dN_standardized, y = omega)) +
  geom_point() +
  xlab(expression("Standardized " * italic(d)[N])) +
  ylab(expression(italic(d)[N]/italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme_set(theme_classic(base_size = 14)); p4
cor.test(x = M0.omega$dN_standardized, y = M0.omega$omega, method = "pearson")

p5 <- ggplot(data = M0.omega, aes(x = dS_standardized, y = omega)) +
  geom_point() +
  xlab(expression("Standardized " * italic(d)[S])) +
  ylab(expression(italic(d)[N]/italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme_set(theme_classic(base_size = 14)); p5
cor.test(x = M0.omega$dS_standardized, y = M0.omega$omega, method = "pearson")

p <- p4/p5
ggsave("~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/dNdSomega.png", plot = p, width = 8, height = 12)


#only core proteome
M0.core.protein <- M0.omega %>%
  filter(number_family_with_loss == 0)

p1 <- ggplot(data = M0.core.protein, aes(y = mean_GC_in_orthogroup, x = dN_standardized)) +
  geom_point() +
  ylab("GC-content") +
  xlab(expression("Standardized " * italic(d)[N])) +  
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.core.protein)), hjust = 1.5, vjust = 2) +
  theme_set(theme_classic(base_size = 14)); p1

p2 <- ggplot(data = M0.core.protein, aes(y = mean_GC_in_orthogroup, x = dS_standardized)) +
  geom_point() +
  ylab("GC-content") +
  xlab(expression("Standardized " * italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.core.protein)), hjust = 1.5, vjust = 2) +
  theme_set(theme_classic(base_size = 14)); p2

p3 <- ggplot(data = M0.core.protein, aes(y = mean_GC_in_orthogroup, x = omega)) +
  geom_point() +
  ylab("GC-content") +
  xlab(expression(italic(d)[N]/italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.core.protein)), hjust = 1.5, vjust = 2) +
  theme_set(theme_classic(base_size = 14)); p3
cor.test(x = M0.core.protein$mean_GC_in_orthogroup, y = M0.core.protein$omega, method = "pearson")


p <- p1/p2/p3
ggsave("~/Desktop/Thesis/Typst/figures/chap3results/GeneByGene/GC_correlations_core_protein.png", plot = p, width = 8, height = 12)


