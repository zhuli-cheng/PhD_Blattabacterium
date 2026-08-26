source("Main.R")

#remove the gene (amiD) with omega > 1
M0.omega <- M0.omega %>%
  filter(omega < 1)


#omega
p1 <- ggplot(data = M0.omega, aes(x = dN_standardized)) +
  theme_classic() +
  geom_density() + 
  xlab(expression("Standardized " * italic(d)[N])) +
  ylab("Density") +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme(text = element_text(size = 14)); p1

mean(M0.omega$dN_standardized)
sd(M0.omega$dN_standardized)
range(M0.omega$dN_standardized)
sum(M0.omega$dN_standardized < 1)

p2 <- ggplot(data = M0.omega, aes(x = mean_position_in_orthogroup, y = dN_standardized)) +
  theme_classic() +
  geom_point() + 
  xlab("Chromosomal position") +
  ylab(expression("Standardized " * italic(d)[N])) +
  annotate("text", x = Inf, y = Inf, label = paste0("N = ", nrow(M0.omega)), hjust = 1.5, vjust = 2) +
  theme(text = element_text(size = 14)); p2
cor.test(M0.omega$dN_standardized, M0.omega$mean_position_in_orthogroup, method = "spearman")

