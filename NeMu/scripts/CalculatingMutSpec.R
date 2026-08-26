source("Main.R")

mutation.spectrum.ff.list <- list()
ms12ff.list <- list()
cosine.similarity <- data.frame(
  group = character(),
  cosine_similarity = numeric(),
  stringsAsFactors = FALSE
)

for (group in selected.groups$group) {
  observed.mutations.192 <- read.table(paste0("../input/monophyleticGroups/NeMuPipeline/", group, "/output/tables/observed_mutations.tsv"), header = T) 
  
  #only ff
  observed.mutations.192.ff <- observed.mutations.192 %>%
    filter(Label == "2")
  observed.mutations.12.ff <- observed.mutations.192.ff %>%
    extract(Mut, into = c("OriNt", "DestNt"), regex = ".*\\[([A-Z])>([A-Z])\\].*") %>%
    mutate(Mut = paste0(OriNt, ">", DestNt))

  mutation.spectrum.ff <- observed.mutations.12.ff %>%
    select(OriNt, DestNt, ProbaFull, Mut) %>%
    group_by(Mut) %>%
    summarise(NumberSubstitution=sum(ProbaFull), .groups = "drop")

  #load program results for comparison
  ms12ff <- read.table(paste0("../input/monophyleticGroups/NeMuPipeline/", group, "/output/tables/ms12ff.tsv"), header = TRUE)
 
  cosine <- coop::cosine(mutation.spectrum.ff$NumberSubstitution, ms12ff$ObsNum)
  
  mutation.spectrum.ff.list[[group]] <- mutation.spectrum.ff
  ms12ff.list[[group]] <- ms12ff
  cosine.similarity <- rbind(cosine.similarity, data.frame(
    group = group,
    cosine_similarty = cosine))

  #plot observed and expected substitutions from NeMu
  NeMu.order <- c(4, 9, 5, 8, 6, 7, 10, 3, 11, 2, 12, 1)
  ms12ff.ordered <- ms12ff[NeMu.order, ]
  ms12ff.ordered$Mut <- factor(
    ms12ff.ordered$Mut,
    levels = ms12ff.ordered$Mut)
  
  p.expected <- ggplot(ms12ff.ordered, aes(x = Mut, y = ExpNum)) +
    geom_col() +
    labs(x = NULL, y = NULL) +
    theme(axis.text  = element_text(size = 14, face = "bold"))
  ggsave(paste0("~/Desktop/Thesis/Typst/figures/chap3results/MutationSpectrum/ExpectedSubstitutions/", group, ".png"), plot = p.expected, width = 6, height = 4, dpi = 300)
  
  p.observed <- ggplot(ms12ff.ordered, aes(x = Mut, y = ObsNum)) +
    geom_col() +
    labs(x = NULL, y = NULL) +
    theme(axis.text  = element_text(size = 14, face = "bold"))
  ggsave(paste0("~/Desktop/Thesis/Typst/figures/chap3results/MutationSpectrum/ObservedSubstitutions/", group, ".png"), plot = p.observed, width = 6, height = 4, dpi = 300)
}
#the calculated spectrum correspond to the NeMu output based on cosine similarity (>99%)

#calculate cosine similarity between any pair of two groups of samples
mat <- sapply(ms12ff.list, function(df) df$MutSpec)
cosine.matrix <- as.data.frame(coop::cosine(mat)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

cosine.matrix <- as.matrix(cosine.matrix)
cosine.matrix[upper.tri(cosine.matrix)] <- ""
cosine.matrix <- as.data.frame(cosine.matrix)



write.table(
  cosine.matrix,
  "~/Desktop/Thesis/Typst/figures/chap3results/MutationSpectrum/spectra_cosine_similarity.csv",
  sep = ",",
  row.names = TRUE,
  col.names = FALSE,
  quote = FALSE)
