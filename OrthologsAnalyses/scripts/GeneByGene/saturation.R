orthogroup <- "573"

# read alignment
dna <- read.dna(paste0("../244/input/alignments/orthogroup_", orthogroup, "_trimmed.ffn"), format = "fasta")

# function to extract codon positions
get_pos <- function(x, pos) {
  x[, seq(pos, ncol(x), by = 3)]
}

dna_pos1 <- get_pos(dna, 1)
dna_pos2 <- get_pos(dna, 2)
dna_pos3 <- get_pos(dna, 3)

# function to compute values
compute_tstv <- function(dna) {
  list(
    dist = dist.dna(dna, model = "raw"),
    ts   = dist.dna(dna, model = "TS"),
    tv   = dist.dna(dna, model = "TV")
  )
}

res1 <- compute_tstv(dna_pos1)
res2 <- compute_tstv(dna_pos2)
res3 <- compute_tstv(dna_pos3)

plot_tstv <- function(res, title) {
  plot(as.vector(res$dist), as.vector(res$ts),
       pch = 16, col = rgb(0,0,1,0.3),
       xlab = "Genetic distance (p-distance)",
       ylab = "Substitutions",
       main = title)
  
  points(as.vector(res$dist), as.vector(res$tv),
         pch = 16, col = rgb(1,0,0,0.3))
  
  legend("topleft", legend = c("Transitions", "Transversions"),
         col = c("blue", "red"), pch = 16)
}

par(mfrow = c(1,3))
plot_tstv(res1, "1st codon position")
plot_tstv(res2, "2nd codon position")
plot_tstv(res3, "3rd codon position")


