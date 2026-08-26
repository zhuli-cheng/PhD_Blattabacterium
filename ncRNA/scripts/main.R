library(readr)
library(dplyr)
library(data.table)
library(readxl)
library(tidyr)

#copied input from OESannot _predict infernal.out files

#define a function to read the infernal ncRNA files
infernal_read <- function(file) {
  df <- fread(file, skip = 2)
  colnames(df) <- c(
    "idx",
    "target_name",
    "accession",
    "query_name",
    "accession1",
    "clan_name",
    "mdl",
    "mdl_from",
    "mdl_to",
    "sequence_from",
    "sequence_to",
    "strand",
    "trunc",
    "pass",
    "gc",
    "bias",
    "score",
    "E_value",
    "inc",
    "olp",
    "anyidx",
    "afrct1",
    "afrct2",
    "winidx",
    "wfrct1",
    "wfrct2",
    "description_of_target"
  )
  df$sample <- tools::file_path_sans_ext(basename(file))
  df$sample <- sub("_infernal.*$", "", df$sample)
  return(df)
}

#list the 316 files
files <- list.files(
  "../input",
  pattern = "infernal\\.out$",
  full.names = TRUE
)


ncRNA_all <- do.call(rbind, lapply(files, infernal_read))

#load sample family
taxonomy <- read.csv("~/Desktop/BlattabacteriumRScripts/TaxonTips/input/host.csv") %>%
  mutate(family_combined = if_else (host.family != "" & !is.na(host.family), host.family, host.family.inferred.from.phylogeny))

ncRNA_all <- ncRNA_all %>%
  left_join(taxonomy, by = join_by(sample == Sample)) %>%
  select(!c(host.family, host.family.inferred.from.phylogeny, Species.determination))

ncRNA_all <- ncRNA_all %>%
  select(target_name, accession, mdl, mdl_from, mdl_to, sequence_from, sequence_to, strand, trunc, gc, score, E_value, sample, family_combined)
