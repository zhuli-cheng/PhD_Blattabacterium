source("LoadPackages.R")

groups <- read_excel("../input/monophyleticGroups/groups/groups.xlsx", sheet = "Sheet1")

ff.count <- ""
for (group in groups$group) {
  observed.mutations.192 <- read.table(paste0("../input/monophyleticGroups/NeMuPipeline/", group, "/output/tables/observed_mutations.tsv"), header = T) 
  
  observed.mutations.192.ff <- observed.mutations.192 %>%
    filter(Label == "2")
  
  ff.count[group] <- nrow(observed.mutations.192.ff)
}
