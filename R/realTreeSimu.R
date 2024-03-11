# TODO with real tree
# Source : https://genome.cshlp.org/content/29/1/53
# Normaliser l'arbre
require("ape")
require("phylolm")
require("phytools")
require("here")
library(tidyverse)
library(ggplot2)
library(patchwork)

source(here("R","utils.R"))

K <- 2

nb_species <- 43

plot_group_on_tree <- function(tree, groups) {
    plot(tree)
    tiplabels(bg = groups, pch = 21)
}

tree <- read.tree(here("R","chen2019.tree"))
# Normalising tree edge length
taille_tree <- diag(vcv(tree))[1]
tree$edge.length <- tree$edge.length / taille_tree

plotTree(tree, ftype="i", node.numbers = TRUE)

# Mus et Rat vs le reste

group_mus_rat_vs_other <- sapply(1:nb_species, function(tip) {
    if (tip %in% getDescendants(tree = tree, 55)) {
        return(1)
    }
    return(2)
})

random_groups <- sample(1:K, nb_species, replace = TRUE)

plot_group_on_tree(tree, group = group_mus_rat_vs_other)
plot_group_on_tree(tree, group = random_groups)


# Generate data for rat&mus vs the rest
N <- 500
base_values <- c(0, 1)
sigma2_phylo <- 1
sigma2_measure <- 0.1
risk_threshold <- 0.05

## Standardized parameters
total_variance <- 1.0 # sigma2_phylo + sigma2_error, fixed [as tree_height = 1]
heri <- c(0.3, 0.5, 0.7, 0.9) # heritability her = sigma2_phylo / total_variance. 0 means only noise. 1 means only phylo.
snr <- 1 # signal to noise ratio snr = size_effect / total_variance

ggsave <- function(..., bg = "white") ggplot2::ggsave(..., bg = bg)
for (her in heri) {
    sim <- N_simulation_typeI_power(N,
        groups_list = list(ratmus_vs_other =  group_mus_rat_vs_other, 
        random = random_groups),
        base_values = c(0, snr * total_variance), 
        sigma2_phylo = her * total_variance,
        sigma2_measure = (1 - her) * total_variance,
        REML = TRUE
    )



    df_sim_plot <- compute_power_typeI(df = sim)

    res_sim_plot <- plot_method_comparison(df_sim_plot, title = paste("Rat&Mus héritabilité ", her))
    res_sim_plot
    ggsave(paste0("img/ratmus_vs_other_her_", her, ".png"), plot = res_sim_plot)
}