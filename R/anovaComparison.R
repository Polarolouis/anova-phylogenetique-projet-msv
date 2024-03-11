# Phylocomparison tools
library(phylolm)
library(phylotools)
library(phytools)
library(phylolimma)
library(ape)
library(tidyverse)

# Plotting
library(ggplot2)
library(patchwork)

# Sourcing the utils
source("./R/utils.R")

# Fixing randomness for reproducibility
set.seed(1234)

# Parameters
nb_species <- 20

# Generating the phylo tree
tree <- rphylo(nb_species, birth = 0.1, death = 0)

# Group selections tries to have group of same size
plotTree(tree, node.numbers = TRUE)

# Here I chose two ancestors to split in two the tree
ancestors <- c(22, 23)
K <- length(ancestors) # The number of groups

# I assign the groups numbers
## Matching the phylogeny
phylomatching_groups <- sapply(1:nb_species, function(tip) {
    get_phylo_group(tip,
        tree,
        ancestors = ancestors
    )
})

## Randomly
random_groups <- sample(1:K, nb_species, replace = TRUE)

## Randomly but with same size of groups
# sameSize_random_groups <- sample(1:K,
#     nb_species,
#     replace = TRUE,
#     prob = table(phylo_matching_groups)
# )
# group_sizes <- table(phylo_matching_groups)


# Saving images of tree
plot_group_on_tree <- function(tree, groups) {
    plot(tree, show.tip.label = FALSE)
    tiplabels(bg = groups, pch = 21)
}

# Saving trees
png(file = "img/group_phylo_matching_tree.png")
plot_group_on_tree(tree, group = phylomatching_groups)
dev.off()

png(file = "img/group_random_tree.png")
plot_group_on_tree(tree, group = random_groups)
dev.off()

# Normalising tree edge length
taille_tree <- diag(vcv(tree))[1]
tree$edge.length <- tree$edge.length / taille_tree

N <- 500
base_values <- c(0, 1)
sigma2_phylo <- 1
sigma2_measure <- 0.1
risk_threshold <- 0.05

## Standardized parameters
total_variance <- 1.0 # sigma2_phylo + sigma2_error, fixed [as tree_height = 1]
heri <- c(0.0, 0.5, 0.75, 1.0) # heritability her = sigma2_phylo / total_variance. 0 means only noise. 1 means only phylo.
snr <- 1 # signal to noise ratio snr = size_effect / total_variance

## Try several parameter values
ggsave <- function(..., bg = "white") ggplot2::ggsave(..., width = 6, height = 6,bg = bg)
for (her in heri) {
    groups_list <- list(phylo = phylomatching_groups, 
        random = random_groups)
    sim <- N_simulation_typeI_power(N,
        groups_list = groups_list,
        base_values = c(0, snr * total_variance), 
        sigma2_phylo = her * total_variance,
        sigma2_measure = (1 - her) * total_variance,
    )



    df_sim_plot <- compute_power_typeI(df = sim)

    res_sim_plot <- plot_method_comparison(df_sim_plot, title = paste("BM héritabilité ", her))
    res_sim_plot
    ggsave(paste0("img/simulation_BM_her_", her, ".png"), plot = res_sim_plot)
}
