source("Main.R")
library(Biostrings)


########## ancestral GC
group <- "group3"

expected.frequencies <-  read.table(paste0("../input/monophyleticGroups/NeMuPipeline/", group, "/output/tables/expected_freqs.tsv"), header = T) 

#MRCA <- "Node9"

ancestral.GC <- expected.frequencies %>%
#  filter(Node == MRCA) %>%
  filter(Label == "all") %>%
  select(Node, Label, A.C, C.A, G.A, T.A) %>%
  mutate(all_nucleotides = sum(A.C, C.A, G.A, T.A)) %>%
  mutate(GC = sum(C.A, G.A)) %>%
  mutate(GC_content = GC/all_nucleotides)

######### sample GC
alignment <- readDNAStringSet("../input/monophyleticGroups/alignments/group3_FcC_supermatrix.fas")
gc_counts <- letterFrequency(alignment, letters = c("G", "C"))
total_bases <- letterFrequency(alignment, letters = c("A", "T", "G", "C"))

gc_no_gaps <- rowSums(gc_counts) / rowSums(total_bases)

gc_no_gaps
