source("LoadPackages.R")

#load groups
groups <- read_excel("../input/monophyleticGroups/groups/groups.xlsx", sheet = "Sheet1")

groups.typst <- groups %>%
  select(group, host_family, number_of_samples, outgroup)

write.table(
  groups.typst,
  "~/Desktop/Thesis/Typst/figures/chap2methods/NeMu_groups.csv",
  sep = ",",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE)


#list of groups selected
selected.groups <- groups %>%
  filter(confidence == "yes") 

selected.groups.typst <- selected.groups %>%
  select(-confidence) %>%
  select(-outgroup) %>%
  select(-rDNAs_location)

#write.table(
#  selected.groups.typst,
#  "~/Desktop/Thesis/Typst/figures/chap3results/MutationSpectrum/groups.csv",
#  sep = ",",
#  row.names = FALSE,
#  col.names = FALSE,
#  quote = FALSE)






