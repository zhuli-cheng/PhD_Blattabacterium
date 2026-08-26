source("Main.R")

spectrum <- list()
names.vec <- c()
k <- 1

for (group in selected.groups$group) {
  types <- c("ff")
  
  for (type in types) {
    ms12 <- read.table(paste0("../input/monophyleticGroups/NeMuPipeline/", group, "/output/tables/ms12", type, ".tsv"), header = TRUE)
    
    ms12 <- ms12 %>%
      select(Mut, ObsNum, ExpNum, MutSpec) %>%
      separate(Mut, into = c("OriginalNucleotide", "DestinationNucleotide"), sep = ">")
    
    M <- matrix(0, 4, 4, dimnames = list(c("A","C","G","T"), c("A","C","G","T")))
    
    for (i in seq_len(nrow(ms12))) {
      from <- ms12$OriginalNucleotide[i]
      to   <- ms12$DestinationNucleotide[i]
      M[from, to] <- ms12$MutSpec[i]}
    
    diag(M) <- -rowSums(M)
    
    vec <- M[row(M) != col(M)]
    
    spectrum[[k]] <- vec
    names.vec[k] <- group
    k <- k + 1
  }
}

names(spectrum) <- names.vec
spec.matrix <- do.call(rbind, spectrum)

cos.mat <- as.matrix(simil(spec.matrix, method = "cosine"))
