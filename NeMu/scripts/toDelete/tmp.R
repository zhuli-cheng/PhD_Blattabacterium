tmp <- frequencies.all %>%
  filter(mutation.type == "nonsyn") %>%
  select(empirical.GC, GC.equilibrium) %>%
  unique()

plot(tmp$empirical.GC ~ tmp$GC.equilibrium)
abline(0, 1)
