#run correlations.R and GC123.R

figure2 <- ((p.correlation1 / p.correlation2 / p.correlation3) | (p.correlation4 / p.correlation5 / p.correlation6)) &
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "plain"), plot.tag.position = c(0.15, 0.9), plot.tag.location = "panel"); figure2


ggsave(filename = "../output/figure2.png", plot = figure2, width = 12, height = 12, dpi = 1200)
