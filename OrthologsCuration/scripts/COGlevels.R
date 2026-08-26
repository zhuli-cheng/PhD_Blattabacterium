#find out how each orthogroup corresponds to each eggNOG_OG at OG_root level (Bacteria level)
root.Blattabacteriaceae <- dataset %>%
  count(orthogroup, OG_root, OG_Blattabacteriaceae, Description, name = "no_of_genes")

Flavobacteriia.Blattabacteriaceae <- dataset %>%
  count(orthogroup, OG_Flavobacteriia, OG_Blattabacteriaceae, Description, name = "no_of_genes")

#check if all eggNOG_OG at OG_root correspond to eggNOG_OG at OG_Blattabacteriaceae level
View(root.Blattabacteriaceae %>%
       group_by(OG_root) %>%
       mutate(unique_Bs = n_distinct(OG_Blattabacteriaceae)) %>%
       ungroup() %>%
       group_by(OG_Blattabacteriaceae) %>%
       mutate(unique_OG_roots = n_distinct(OG_root)) %>%
       ungroup())

View(Flavobacteriia.Blattabacteriaceae %>%
       group_by(OG_Flavobacteriia) %>%
       mutate(unique_Bs = n_distinct(OG_Blattabacteriaceae)) %>%
       ungroup() %>%
       group_by(OG_Blattabacteriaceae) %>%
       mutate(unique_OG_Flavobacteriia = n_distinct(OG_Flavobacteriia)) %>%
       ungroup())
#because the same OG at the B level always give 1 or at most 2 OGs at the OG_root level. And when it is 2, the functions are still the same. Therefore we will use OGs at B level. 
