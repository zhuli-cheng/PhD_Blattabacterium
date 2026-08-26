source("Main.R")

frequency <- list()
k <- 1

for (group in selected.groups$group) {
  types <- c("nonsyn", "syn", "ff")
  
  for (type in types) {
    #read input
    ms12 <- read.table(paste0("../input/monophyleticGroups/NeMuPipeline/", group, "/output/tables/ms12", type, ".tsv"), header = T) 
    ms12 <- ms12 %>%
      select(Mut, ObsNum, ExpNum, MutSpec) %>%
      separate(Mut, into = c("OriginalNucleotide", "DestinationNucleotide"), sep = ">")
    
    #create mutation spectrum matrix and enter data
    M <- matrix(0, 4, 4,
                dimnames = list(c("A","C","G","T"),
                                c("A","C","G","T")))
    
    for (i in seq_len(nrow(ms12))) {
      from <- ms12$OriginalNucleotide[i]
      to   <- ms12$DestinationNucleotide[i]
      M[from, to] <- ms12$MutSpec[i]
    }
    
    diag(M) <- -rowSums(M)
    
    #calculate stationary distribution
    e <- eigen(t(M))
    v <- e$vectors
    w <- e$values
    
    stat <- Re(v[, which.min(abs(w))])
    pi <- stat / sum(stat)
    
    #empirical frequencies
    empirical <- ms12 %>%
      group_by(OriginalNucleotide) %>%
      mutate(substitution.sum = sum(ObsNum)) %>%
      mutate(site.sum = sum(ExpNum)) %>%
      ungroup() %>%
      select(OriginalNucleotide, substitution.sum, site.sum) %>%
      unique() %>%
      mutate(frequency = site.sum/sum(site.sum)) %>%
      mutate(empirical.GC = frequency[OriginalNucleotide == "G"] + frequency[OriginalNucleotide == "C"])
      
    #global GC content per genome
     global <- read.table(paste0("../input/monophyleticGroups/NeMuPipeline/", group, "/output/tables/ms12all.tsv"), header = T) %>%
      separate(Mut, into = c("OriginalNucleotide", "DestinationNucleotide"), sep = ">") %>%
      select(OriginalNucleotide, ExpNum) %>%
      unique() %>%
      mutate(frequency = ExpNum/sum(ExpNum)) %>%
      mutate(global.GC = frequency[OriginalNucleotide == "G"] + frequency[OriginalNucleotide == "C"])

    #combine results into a data frame
    frequencies <- data.frame(
      nucleotide = c("A", "C", "G", "T"),
      group = group,
      family = groups$host_family[groups$group == group],
      mutation.type = type,
      substitutions = empirical$substitution.sum,
      stationary.distribution = pi,
      empirical.frequency = empirical$frequency,
      GC.equilibrium = as.numeric(pi[2] + pi[3]),
      empirical.GC = empirical$empirical.GC,
      global.GC = global$global.GC
    )
    
    frequency[[k]] <- frequencies
    k <- k + 1
  }
}

frequencies.all <- do.call(rbind, frequency)

write.table(
  frequencies.all,
  "../output/monophyleticGroups/StationaryDistribution.csv",
  sep = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE)



##############
#double-check whether the frequencies are stationary (frequency of change sum = 0)
names(pi) <- c("A", "C", "G", "T")
#A
pi["A"]*M["A", "C"] + pi["A"]*M["A", "G"] + pi["A"]*M["A", "T"] - pi["C"]*M["C", "A"] - pi["G"]*M["G", "A"] - pi["T"]*M["T", "A"]
#C
pi["C"]*M["C", "A"] + pi["C"]*M["C", "G"] + pi["C"]*M["C", "T"] - pi["A"]*M["A", "C"] - pi["G"]*M["G", "C"] - pi["T"]*M["T", "C"]
#G
pi["G"]*M["G", "A"] + pi["G"]*M["G", "C"] + pi["G"]*M["G", "T"] - pi["A"]*M["A", "G"] - pi["C"]*M["C", "G"] - pi["T"]*M["T", "G"]
#T
pi["T"]*M["T", "A"] + pi["T"]*M["T", "C"] + pi["T"]*M["T", "G"] - pi["A"]*M["A", "T"] - pi["C"]*M["C", "T"] - pi["G"]*M["G", "T"]

