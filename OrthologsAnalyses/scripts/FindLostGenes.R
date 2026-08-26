group1 <- sample.orthogroup %>%
  filter(sample == "CAM19_2_C136_CR_S38_L004" | sample == "CAM19_1_C025_CR_S25_L001" | sample == "CAM19_2_C108_CR_S11_L004") 

rownames(group1) <- group1$sample

group1 <- group1[, group1["CAM19_2_C136_CR_S38_L004", ] != 0 & group1["CAM19_1_C025_CR_S25_L001", ] == 0 & group1["CAM19_2_C108_CR_S11_L004", ] == 0]

group1.orthogroups <- colnames(group1)
group1.orthogroups <- gsub("new_orthogroup", "", group1.orthogroups)

group1.orthogroups <- orthogroup.sample %>%
  filter(new_orthogroup %in% group1.orthogroups)



group2 <- sample.orthogroup %>%
  filter(sample == "CAM19_2_C149_CR_S49_L004" | sample == "AU076_S405") 

rownames(group2) <- group2$sample

group2 <- group2[, group2["CAM19_2_C149_CR_S49_L004", ] == 0 & group2["AU076_S405", ] != 0 ]

group2.orthogroups <- colnames(group2)
group2.orthogroups <- gsub("new_orthogroup", "", group2.orthogroups)

group2.orthogroups <- orthogroup.sample %>%
  filter(new_orthogroup %in% group2.orthogroups)




group3 <- sample.orthogroup %>%
  filter(sample == "CAM19_2_C149_CR_S49_L004" | sample == "AU076_S405") 

rownames(group3) <- group3$sample

group3 <- group3[, group3["CAM19_2_C149_CR_S49_L004", ] != 0 & group3["AU076_S405", ] == 0 ]

group3.orthogroups <- colnames(group3)
group3.orthogroups <- gsub("new_orthogroup", "", group3.orthogroups)

group3.orthogroups <- orthogroup.sample %>%
  filter(new_orthogroup %in% group3.orthogroups)

