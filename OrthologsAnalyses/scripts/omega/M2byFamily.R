source("Main.R") 
source("CollapsingHighlySimilarGenomes.R")
collapsed.tree$edge.length[] <- 0

#find families that have gene loss
family.omega <- data.frame(tip_label = collapsed.tree$tip.label) %>%
  mutate(host_family = taxonomy$family[match(collapsed.tree$tip.label, taxonomy$sample)]) 

families <- unique(family.omega$host_family)

rename_branches <- function(families, index) {
  family <- families[index]
  tips <- family.omega$tip_label[family.omega$host_family == family]
  ancestral.node <- getMRCA(collapsed.tree, tips)
  descendant_edges <- which(collapsed.tree$edge[, 1] == ancestral.node | collapsed.tree$edge[, 2] %in% getDescendants(collapsed.tree, ancestral.node, "all"))
  collapsed.tree$edge.length[descendant_edges] <<- index
}

for (i in 1:length(families)) {
  rename_branches(families, i)
}

collapsed.tree$node.label <- NULL
ggtree(collapsed.tree) %<+% family.omega +
  geom_tiplab(size = 2, hjust = -0.1, aes(label = host_family))

write.tree(collapsed.tree, "../244/output/M2_100.tree")

########### run M2 models on Deigo, then: #############

################ analysis
dS.tree <- read.tree("tmp/dS.tree")
plot(dS.tree)
ggtree(dS.tree) %<+% taxonomy.tip +
  geom_tiplab(size = 2, hjust = -0.1, aes(label = taxon))

dN.tree <- read.tree("tmp/dN.tree")
plot(dN.tree)
ggtree(dN.tree) %<+% taxonomy.tip +
  geom_tiplab(size = 2, hjust = -0.1, aes(label = taxon))

w.tree <- read.tree("tmp/w.tree")
plot(w.tree)


#load eggNOG annotations
eggNOG <- read.csv("../../OrthologsCuration/244/input/383_diamond.emapper.annotations",
                   sep = "\t",
                   comment.char = "#",
                   header = F,
                   na.strings = "-")
colnames(eggNOG) <- c("query", "seed_ortholog", "evalue", "score", "eggNOG_OGs",
                      "max_annot_lvl", "COG_category", "Description", "Preferred_name",
                      "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction",
                      "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy",
                      "BiGG_Reaction", "PFAMs")


#create table showing the loss of genes in each families
combined_family_presence <- data.frame()
for (i in 1:nrow(orthogroup.sample)) {
  #load gene presence/absence
  orthogroup.name <- paste0("new_orthogroup", i)
  orthogroup.presence <- setNames(ifelse(sample.orthogroup[[orthogroup.name]] == 0, 0, 1), sample.orthogroup$sample)
  
  #find families that have gene loss
  family.presence <- data.frame(tip_label = collapsed.tree$tip.label) %>%
    mutate(family = taxonomy$family[match(collapsed.tree$tip.label, taxonomy$sample)]) %>%
    mutate(gene_presence = orthogroup.presence[collapsed.tree$tip.label]) %>%
    group_by(family) %>%
    mutate(family_loss = ifelse(sum(gene_presence) == 0, "complete_gene_loss", ifelse(sum(gene_presence) == n(), "no_gene_loss", "partial_gene_loss"))) %>%
    ungroup() %>%
    select(family, family_loss) %>%
    unique() %>%
    mutate(new_orthogroup = i)
  
  combined_family_presence <- rbind(combined_family_presence, family.presence)
}

#load omega data
M2.omega <- read.csv("../244/input/M2_summary.txt") %>%
  pivot_longer(
    cols = -1,
    names_to = "family",
    values_to = "omega"
  ) 

M2.omega <- M2.omega %>%
  mutate(family = ifelse(family == "w1", "Blaberidae", family)) %>%
  mutate(family = ifelse(family == "w2", "Blattellidae", family)) %>%
  mutate(family = ifelse(family == "w3", "Anaplectidae", family)) %>%
  mutate(family = ifelse(family == "w4", "Blattidae", family)) %>%
  mutate(family = ifelse(family == "w5", "Corydiidae", family)) %>%
  mutate(family = ifelse(family == "w6", "Pseudophyllodromiidae", family)) %>%
  mutate(family = ifelse(family == "w7", "Ectobiidae", family)) 

M2.omega <- M2.omega %>%
  left_join(orthogroup.sample, by = "new_orthogroup") %>%
  select(new_orthogroup, family, omega, curated_COG_category, OGs_in_orthogroup, gene, mean_GC_in_orthogroup, mean_length_in_orthogroup)

M2.omega <- full_join(combined_family_presence, M2.omega, by = c("new_orthogroup", "family"))

M2.omega <- M2.omega %>%
  mutate(across(c("family", "family_loss", "new_orthogroup", "curated_COG_category"), factor))

ggplot(M2.omega, aes(x = family, y = omega)) +
  geom_boxplot(fill = "skyblue", color = "darkblue") +
  theme_minimal() +
  labs(
    x = "Host family",
    y = "Omega"
  )

ggplot(M2.omega, aes(x = curated_COG_category, y = omega)) +
  geom_boxplot(fill = "skyblue", color = "darkblue") +
  theme_minimal() +
  labs(
    x = "COG categories",
    y = "Omega"
  )

ggplot(M2.omega, aes(x = family_loss, y = omega)) +
  geom_boxplot(fill = "skyblue", color = "darkblue") +
  theme_minimal() +
  labs(
    x = "Gene loss within each family",
    y = "Omega"
  )


#GLMM
M2.omega.binary <- M2.omega %>%
  filter(family_loss != "complete_gene_loss") %>%
  filter(!is.na(omega)) %>%
  filter(mean_length_in_orthogroup > 600)

M2.omega.binary <- M2.omega.binary %>%
  group_by(new_orthogroup) %>%
  mutate(mean_omega_in_orthogroup = mean(omega)) %>%
  mutate(sd_omega_in_orthogroup = sd(omega)) %>%
  mutate(deviation_omega_in_orthogroup = (omega - mean(omega))/sd(omega)) %>%
  ungroup() %>%
  filter(!is.na(deviation_omega_in_orthogroup))

model <- glmer(family_loss ~ 
               omega +
                 (1 | new_orthogroup) +
                 (1 | family),
               data = M2.omega.binary,
               family = binomial)

#assumptions
simout  <-  simulateResiduals (model, n = 250); plot(simout)
hist(intercep.nb <- coef(model)$new_orthogroup[, 1], breaks = 100)
hist(intercep.nb <- coef(model)$family[, 1], breaks = 100)
hist(model@u,  breaks = 100)

#resumodel#results
summary(model)
#model result
fixed <- fixef(model); fixed
confintfixed <- confint(model, parm = "beta_", method = "Wald"); confintfixed

#plot
M2.omega.binary$predicted <- predict(model, type = "response")
M2.omega.binary$raw <- ifelse(M2.omega.binary$family_loss == "no_gene_loss", 0, 1)
ggplot(M2.omega.binary, aes(x = omega)) +
  geom_point(aes(y = raw), alpha = 0.6, color = "blue") +
  stat_smooth(inherit.aes = F, aes(x = omega, y = predicted),
              method = "glm",
              method.args = list(family = "binomial"), #replace with whatever family you use
              color = "black", size = 0.5, level = 0.95,
              se = TRUE) +  
  theme_minimal() +
  labs(x = "omega", y = "gene loss within family")

M2.omega <- M2.omega %>%
  group_by(new_orthogroup) %>%
  filter(!is.na(omega)) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  filter(n == 7)


model <- lm(data = M2.omega, omega ~ as.factor(family) + as.factor(new_orthogroup)) 
summary(model)
