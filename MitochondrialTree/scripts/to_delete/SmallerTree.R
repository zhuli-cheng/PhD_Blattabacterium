source("Main.R")

#full tree
ID.tip.tree <- ggtree(mito.tree, ladderize = FALSE) %<+% absence +
  geom_tiplab(size = 2.5, hjust = -0.1, aes(color = absence_confirmed)) +  
  geom_treescale(x = 0, y = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +
  geom_text2(aes(subset = !isTip, label = label),
             hjust = -0.3,
             vjust = -0.3,
             size = 2); ID.tip.tree
#ggsave(filename = "../output/BlattabacteriumLossTree.pdf", plot = ID.tip.tree, width = 20, height = 50, dpi = 1200, limitsize = FALSE)



modify_tip(mito.tree, "mixed")

absence <- absence %>% 
  full_join(taxonomy.tip, join_by(sample == tip_label)) %>%
  mutate(taxon = ifelse(is.na(taxon) | taxon == "", sample, taxon))

absence.mixed <- absence %>%
  mutate(absence_confirmed = ifelse(is.na(absence_confirmed), 0, absence_confirmed)) %>%
  mutate(absence_confirmed = as.factor(absence_confirmed))

absence.mixed$taxon[absence.mixed$sample == "CAM19_2_C164_CR_S63_L004"] <- paste0(absence.mixed$taxon[absence.mixed$sample == "CAM19_2_C164_CR_S63_L004"], " (contamination?)")
absence.mixed$taxon[absence.mixed$sample == "CAM19_1_C037_CR_S37_L001"] <- paste0(absence.mixed$taxon[absence.mixed$sample == "CAM19_1_C037_CR_S37_L001"], " (contamination?)")
absence.mixed$taxon[absence.mixed$sample == "CAM19_2_C007_CR_S9_L002"] <- paste0(absence.mixed$taxon[absence.mixed$sample == "CAM19_2_C007_CR_S9_L002"], " (contamination?)")


ID.tip.tree <- ggtree(mito.tree, ladderize = FALSE) %<+% absence.mixed +
  geom_tiplab(size = 2.5, hjust = -0.1, aes(label = taxon, color = absence_confirmed)) +
  geom_treescale(x = 0, y = 0) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.2))) +
  geom_text2(aes(subset = !isTip, label = label),
             hjust = -0.3,
             vjust = 0,
             size = 2) +
  scale_color_manual(values = c("0" = "black", "1"  = "blue")) + 
  theme(legend.position = "none"); ID.tip.tree
ggsave(filename = "../output/BlattabacteriumLossTreeSmall.pdf", plot = ID.tip.tree, width = 10, height = 16, dpi = 1200)
#ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/BlattabacteriumLossTreeSmall.pdf", plot = ID.tip.tree, width = 10, height = 14, dpi = 1200)



