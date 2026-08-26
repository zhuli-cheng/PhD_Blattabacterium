#after running GC.R, dNdS.R, dS.R

figure1 <- (p.GC1 + p.GC2) / p.GC5 / (p.omega1 + p.omega2) / p.omega5 / (p.dS1 + p.dS2) &
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "plain"), plot.tag.position = c(0.07, 0.95), plot.tag.location = "panel"); figure1

ggsave(filename = "../output/figure1.png", plot = figure1, width = 12, height = 16, dpi = 1200)
