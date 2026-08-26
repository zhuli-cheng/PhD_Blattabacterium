library(genoPlotR)

alignment <- genoPlotR::read_mauve_backbone("~/Desktop/progressiveMauve/77.xmfa.backbone")
plot_gene_map(dna_segs=alignment$dna_segs, comparisons=alignment$comparisons)

