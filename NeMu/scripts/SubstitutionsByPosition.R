source("Main.R")
#the alignment is positioned by: order of the genes based on their position on the chromosome, and using the coding strand sequences


#mutations by gene position
for (group in selected.groups$group) {
  #load reconstructed mutations
  observed.mutations.192 <- read.table(paste0("../input/monophyleticGroups/NeMuPipeline/", group, "/output/tables/observed_mutations.tsv"), header = T)
  #select ff
  observed.mutations.192.ff <- observed.mutations.192 %>%
    filter(Label == 2)
  #select syn
  observed.mutations.192.syn <- observed.mutations.192 %>%
    filter(Label == 1 | Label == 2)
  #select non-syn
  observed.mutations.192.nonsyn <- observed.mutations.192 %>%
    filter(Label == 0)
  
  #location of the rDNA operon
  rDNAs.position <- selected.groups$rDNAs_location[selected.groups$group == group]
  
  #four-fold
  p.ff <- ggplot(observed.mutations.192.ff, aes(x = PosInGene)) +
    geom_histogram(fill = "lightblue", color = "blue", binwidth = 100) +
    labs(x = NULL, y = NULL) +
    geom_vline(xintercept = rDNAs.position, linetype = "dashed", size = 1) +
    theme(axis.text  = element_text(size = 14, face = "bold"))
  #    annotate("text", x = rDNAs.position, y = 30, label = "rDNA operon", vjust = 1.5, angle = 0, size = 5) 
  ggsave(paste0("~/Desktop/Thesis/Typst/figures/chap3results/MutationSpectrum/SubstitutionsByPosition/ff/", group, ".png"), plot = p.ff, width = 12, height = 2, dpi = 300)
  
  #syn
  p.syn <- ggplot(observed.mutations.192.syn, aes(x = PosInGene)) +
    geom_histogram(fill = "lightblue", color = "blue", binwidth = 100) +
    labs(x = NULL, y = NULL) +
    geom_vline(xintercept = rDNAs.position, linetype = "dashed", size = 1) +
    theme(axis.text  = element_text(size = 14, face = "bold"))
  #    annotate("text", x = rDNAs.position, y = 40, label = "rDNA operon", vjust = 1.5, angle = 0, size = 5) 
  ggsave(paste0("~/Desktop/Thesis/Typst/figures/chap3results/MutationSpectrum/SubstitutionsByPosition/syn/", group, ".png"), plot = p.syn, width = 12, height = 2, dpi = 300)
  
  #nonsyn
  p.nonsyn <- ggplot(observed.mutations.192.nonsyn, aes(x = PosInGene)) +
    geom_histogram(fill = "lightblue", color = "blue", binwidth = 100) +
    labs(x = NULL, y = NULL) +
    geom_vline(xintercept = rDNAs.position, linetype = "dashed", size = 1) +
    theme(axis.text  = element_text(size = 14, face = "bold"))
  #    annotate("text", x = rDNAs.position, y = 30, label = "rDNA operon", vjust = 1.5, angle = 0, size = 5) 
  ggsave(paste0("~/Desktop/Thesis/Typst/figures/chap3results/MutationSpectrum/SubstitutionsByPosition/nonsyn/", group, ".png"), plot = p.nonsyn, width = 12, height = 2, dpi = 300)
  
  #all substitutions
  p.all <- ggplot(observed.mutations.192, aes(x = PosInGene)) +
    geom_histogram(fill = "lightblue", color = "blue", binwidth = 100) +
    labs(x = NULL, y = NULL) +
    geom_vline(xintercept = rDNAs.position, linetype = "dashed", size = 1) +
    theme(axis.text  = element_text(size = 14, face = "bold"))
#    annotate("text", x = rDNAs.position, y = 30, label = "rDNA operon", vjust = 1.5, angle = 0, size = 5) 
  ggsave(paste0("~/Desktop/Thesis/Typst/figures/chap3results/MutationSpectrum/SubstitutionsByPosition/all/", group, ".png"), plot = p.all, width = 12, height = 2, dpi = 300)
  
  print(p.ff)
  print(p.syn)
  print(p.nonsyn)
  print(p.all)
}


################# plot substitutions of each mutation type on the position
for (group in selected.groups$group) {
  observed.mutations.192 <- read.table(paste0("../input/monophyleticGroups/NeMuPipeline/", group, "/output/tables/observed_mutations.tsv"), header = T)
  #select synonymous sites
  observed.mutations.192.syn <- observed.mutations.192 %>%
    filter(Label == 2 | Label == 1)
  #transform into 12 types
  observed.mutations.12.syn <- observed.mutations.192.syn %>%
    mutate(mut = gsub(".*\\[(.*)\\].*", "\\1", Mut))

  rDNAs.position <- selected.groups$rDNAs_location[selected.groups$group == group]
  
  mutation.types <- unique(observed.mutations.12.syn$mut)

  #loop over each group and plot
  for (g in mutation.types) {
    p <- ggplot(filter(observed.mutations.12.syn, mut == g), aes(x = PosInGene)) +
      geom_histogram(fill = "lightblue", color = "blue", binwidth = 100) +
      labs(title = NULL, x = NULL, y = NULL) +
      coord_cartesian(ylim = c(0, 20)) +
      geom_vline(xintercept = rDNAs.position, linetype = "dashed", size = 1) +
      theme(axis.text  = element_text(size = 14, face = "bold"))
      # annotate("text", x = rDNAs.position, y = 15, label = "rDNA operon", vjust = 1.5, angle = 0, size = 5) 
    ggsave(paste0("~/Desktop/Thesis/Typst/figures/chap3results/MutationSpectrum/SubstitutionsByPosition/byType/syn/", group, "/", g, ".png"), plot = p, width = 12, height = 2, dpi = 300)
    print(p)
  }
}

for (group in selected.groups$group) {
  observed.mutations.192 <- read.table(paste0("../input/monophyleticGroups/NeMuPipeline/", group, "/output/tables/observed_mutations.tsv"), header = T)
  #select synonymous sites
  observed.mutations.192.ff <- observed.mutations.192 %>%
    filter(Label == 2)
  #transform into 12 types
  observed.mutations.12.ff <- observed.mutations.192.ff %>%
    mutate(mut = gsub(".*\\[(.*)\\].*", "\\1", Mut))
  
  rDNAs.position <- selected.groups$rDNAs_location[selected.groups$group == group]
  
  mutation.types <- unique(observed.mutations.12.ff$mut)
  
  #loop over each group and plot
  for (g in mutation.types) {
    p <- ggplot(filter(observed.mutations.12.ff, mut == g), aes(x = PosInGene)) +
      geom_histogram(fill = "lightblue", color = "blue", binwidth = 100) +
      labs(title = NULL, x = NULL, y = NULL) +
      coord_cartesian(ylim = c(0, 20)) +
      geom_vline(xintercept = rDNAs.position, linetype = "dashed", size = 1) +
      theme(axis.text  = element_text(size = 14, face = "bold"))
    # annotate("text", x = rDNAs.position, y = 15, label = "rDNA operon", vjust = 1.5, angle = 0, size = 5) 
    ggsave(paste0("~/Desktop/Thesis/Typst/figures/chap3results/MutationSpectrum/SubstitutionsByPosition/byType/ff/", group, "/", g, ".png"), plot = p, width = 12, height = 2, dpi = 300)
    print(p)
  }
}

