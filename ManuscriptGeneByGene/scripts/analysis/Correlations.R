source("Main.R")

#remove the gene (amiD) with omega > 1
M0.omega <- M0.omega %>%
  filter(omega < 1)
ko.categories <- ko.categories %>%
  filter(omega < 1)


############## correlations ################
p.correlation1 <- ggplot(data = M0.omega, aes(y = mean_GC_in_orthogroup, x = dN_standardized)) +
  geom_point() +
  ylab("GC content") +
  xlab(expression("Standardized " * italic(d)[N])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme_AnimalBehaviour(); p.correlation1

p.correlation2 <- ggplot(data = M0.omega, aes(y = mean_GC_in_orthogroup, x = dS_standardized)) +
  geom_point() +
  ylab("GC content") +
  xlab(expression("Standardized " * italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme_AnimalBehaviour(); p.correlation2

p.correlation3 <- ggplot(data = M0.omega, aes(y = mean_GC_in_orthogroup, x = omega)) +
  geom_point() +
  ylab("GC content") +
  xlab(expression(italic(d)[N]/italic(d)[S])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme_AnimalBehaviour(); p.correlation3
cor.test(x = M0.omega$mean_GC_in_orthogroup, y = M0.omega$omega, method = "pearson")



