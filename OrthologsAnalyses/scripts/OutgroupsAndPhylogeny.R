source("Main.R")
#select orthogroups (core proteins) to build a phylogeny

omega <- read.csv(paste0("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsAnalyses/", sample.size, "/input/omega_summary.out"), header = T)

phylogeny <- dataset.clean %>%
  filter(core_proteome == TRUE) %>%
  select(new_orthogroup, sample_size_in_orthogroup, mean_length_in_orthogroup, min_length_in_orthogroup, max_length_in_orthogroup) %>%
  unique() %>%
  filter(max_length_in_orthogroup <= 1.3) %>%
  filter(min_length_in_orthogroup >= 0.8)

#outgroups
#do not use ureC when including outgroups
outgroup.blast <- read.table("../244plus9/input/outgroups_blast.out")
colnames(outgroup.blast) <- c("new_orthogroup", "species", "hit", "identity_percentage", "alignment_length", "mismatch", "gap", "evalue", "bitscore")
outgroup.blast <- outgroup.blast %>%
  group_by(new_orthogroup) %>%
  mutate(mean_identity = mean(identity_percentage)) %>%
  mutate(m_identity = identity_percentage/max(identity_percentage)) %>%
  mutate(mean_alignment = mean(alignment_length)) %>%
  mutate(m_alignment = alignment_length/max(alignment_length)) %>%
  mutate(mean_bitscore = mean(bitscore)) %>%
  mutate(m_bitscore = bitscore/max(bitscore)) %>%
  ungroup() %>%
  left_join(orthogroup.sample[,c("new_orthogroup", "gene", "n")], by = c("new_orthogroup"))

outgroup.headers <- read.csv("../244plus9/input/outgroup_headers.csv", header = FALSE)
colnames(outgroup.headers) <- c("species", "hit", "gene_description")

outgroup.blast <- outgroup.blast %>%
  left_join(outgroup.headers, by = c("species", "hit"))

#filter unreliable matches
outgroup.blast <- outgroup.blast %>%
  filter(m_identity > 0.7) %>%
  filter(m_bitscore > 0.9) %>%
  filter(!grepl("ureC", gene)) %>%
  filter(!grepl("ureA", gene)) 
  
#check outgroup hits
tmp <- outgroup.blast %>%
  group_by(new_orthogroup) %>%
  filter(length(unique(gene_description)) > 1)
length(unique(outgroup.blast$new_orthogroup))
outgroup.blast.output <- outgroup.blast %>%
  select(new_orthogroup, species, hit)
write.table(outgroup.blast.output, "../244/output/outgroup_blast.csv", row.names = F, col.names = F, quote = F, sep = ",")

#select orthogroups for phylogenetics
orthogroups.phylogeny <- outgroup.blast %>%
  left_join(omega[, c("new_orthogroup", "omega")], by = join_by(new_orthogroup)) %>%
  filter(new_orthogroup %in% phylogeny$new_orthogroup) %>%
  group_by(new_orthogroup) %>%
  mutate(n.species = n()) %>%
  ungroup() %>%
  filter(n.species > 8) %>%
  filter(omega < 0.1)
length(unique(orthogroups.phylogeny$new_orthogroup))

#output
phylogeny.headers <- dataset.clean %>%
  filter(new_orthogroup %in% orthogroups.phylogeny) %>%
  select(query, new_orthogroup)
write.table(phylogeny.headers, paste0("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsAnalyses/", sample.size, "/output/phylogeny_headers.csv"), row.names = F, col.names = F, quote = F, sep = ",")

orthogroups.phylogeny <- unique(orthogroups.phylogeny$new_orthogroup)
write.table(orthogroups.phylogeny, paste0("/Users/zhulicheng/Desktop/BlattabacteriumRScripts/OrthologsAnalyses/", sample.size, "/output/outgroups_phylogeny.csv"), row.names = F, col.names = F, quote = F, sep = ",")



