aa <- read.csv("../244/input/codon_usage.csv")

p <- ggplot(data = aa, aes(x = expected, y = observed, label = amino_acid)) +
  geom_point(color = "blue") +       # Scatter plot
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +  # x = y line
  geom_text(vjust = -0.5, hjust = 0.5, size = 3) +  # Add labels near points
  theme_minimal(); p
