source("Main.R")
#source("CollapsingHighlySimilarGenomes.R")

################# parsimony method ##################
# function to assign states on internal nodes, based on tip states and the tree
assign_states <- function(collapsed.tree, node.states) {
  num.internal.nodes <- Nnode(collapsed.tree)
  # Traverse the collapsed.tree from the tips to the root
  for (node in (Ntip(collapsed.tree) + 1):(Ntip(collapsed.tree) + num.internal.nodes)) {
    descendants <- getDescendants(collapsed.tree, node)
    descendant.tips <- descendants[descendants <= Ntip(collapsed.tree)]
    if (all(node.states[descendant.tips] == 0)) {
      node.states[node] <- 0  # Assign state 0 if all descendants are 0
    } else {
      node.states[node] <- 1  # Otherwise, assign state 1
    }
  }
  return(node.states)
}

# Function to count transitions based on the node states table
count_changes <- function(tree, node_states) {
  changes <- 0
  # Loop through each edge in the tree
  for (i in 1:nrow(tree$edge)) {
    node_left <- tree$edge[i, 1]
    node_right <- tree$edge[i, 2]
    state_left <- node_states[node_left]
    state_right <- node_states[node_right]
    if (state_left != state_right) {
      changes <- changes + 1
    }
  }
  return(changes)
}

#loop through the orthogroups and save the results in a gene loss table
gene.loss <- data.frame(new_orthogroup = integer(0), gene_loss = integer(0))

for (i in c(1:436, 438:598)) {
  orthogroup_name <- paste0("orthogroup", i)
  orthogroup_data <- setNames(ifelse(sample.orthogroup[[paste0("new_", orthogroup_name)]] == 0, 0, 1), sample.orthogroup$sample)
  orthogroup_data <- orthogroup_data[collapsed.tree$tip.label]
  
  node.states <- rep(NA, Ntip(collapsed.tree) + Nnode(collapsed.tree))
  node.states[1:Ntip(collapsed.tree)] <- orthogroup_data[collapsed.tree$tip.label]
  
  node.states <- assign_states(collapsed.tree, node.states)
  num_changes <- count_changes(collapsed.tree, node.states)
  
  p <- ggtree(collapsed.tree, ladderize = FALSE) +
    geom_tippoint(aes(color = factor(node.states))) +
    geom_nodepoint(aes(color = factor(node.states))) +
    scale_color_manual(values = c("0" = "#F8766D", "1" = "#00BFC4")) +
    labs(title = paste("Orthogroup #", i, ":", orthogroup.sample$gene[orthogroup.sample$new_orthogroup == i])) 
  
# ggsave(filename = paste0("../", sample.size, "/output/geneLossOnTree/orthogroup", i, ".png"), plot = p, width = 10, height = 10)
  gene.loss <- rbind(gene.loss, data.frame(new_orthogroup = i, gene_loss = num_changes))
}

p <- ggplot(data = gene.loss, aes(x = gene_loss)) +
  geom_histogram( binwidth=1, alpha=0.9) +
  theme_classic() +
  xlab("Number of independant loss") +
  ylab("Count of orthogroups") +
  theme(plot.title = element_text(size=15)); p

# ggsave(filename = "~/Desktop/Thesis/Typst/figures/chap3results/gene_loss_histogram.png", plot = p)


################# likelihood method ################
#specify a matrix where 0 cannot transition to 1
#Q <- matrix(c(0, 0, 1, 0), nrow = 2, byrow = T)
#rownames(Q) <- colnames(Q) <- c("0", "1")
#likelihood.result <- ace(orthogroup1, collapsed.tree, model = Q, type = "discrete", marginal = T)#issue caused by the root.edge length being 0
#plot(collapsed.tree)
#nodelabels(pie = likelihood.result$lik.anc, piecol = c("blue", "red"), cex = 0.2)
