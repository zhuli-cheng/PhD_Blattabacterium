source("Main.R")

position <- dataset.clean %>%
  select(new_orthogroup, mean_position_in_orthogroup) %>%
  unique() %>%
  arrange(mean_position_in_orthogroup) %>%
  select(new_orthogroup)

write.table(position, "../244/output/chromosome_position.csv", col.names = FALSE, row.names = FALSE)
