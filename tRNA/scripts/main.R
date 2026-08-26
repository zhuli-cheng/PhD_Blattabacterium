library(readr)
library(dplyr)
library(data.table)
library(readxl)
library(tidyr)

#define a function to read the tRNA scan files
tRNAscan_read <- function(file) {
  df <- fread(file, skip = 3)
  colnames(df) <- c(
    "seqid",
    "trna_id",
    "start",
    "end",
    "aa",
    "anticodon",
    "intron_start",
    "intron_end",
    "score",
    "note"
  )
  df$sample <- tools::file_path_sans_ext(basename(file))
  df$sample <- sub("_tRNASE.*$", "", df$sample)
  return(df)
}

#list the 316 files
files <- list.files(
  "../input",
  pattern = "tRNASE\\.out$",
  full.names = TRUE
)


tRNA_all <- do.call(rbind, lapply(files, tRNAscan_read))

tRNA_all <- tRNA_all %>%
  filter(is.na(note) | note != "pseudo") %>% #remove pseudogene
  filter(aa != "Undet") %>% #remove false prediction 
  mutate(tRNA = paste(aa, anticodon, sep = "_"))


#load sample family
taxonomy <- read.csv("~/Desktop/BlattabacteriumRScripts/TaxonTips/input/host.csv") %>%
  mutate(family_combined = if_else (host.family != "" & !is.na(host.family), host.family, host.family.inferred.from.phylogeny))

tRNA_all <- tRNA_all %>%
  left_join(taxonomy, by = join_by(sample == Sample)) %>%
  select(!c(host.family, host.family.inferred.from.phylogeny, Species.determination))


