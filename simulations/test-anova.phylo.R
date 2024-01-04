library(phylolm)
library(phylotools)
library(phytools)
library(ape)
set.seed(1234)

N <- 100 # Number of different simulations
n <- 100

# Arbre
tree <- rphylo(n, 0.1, 0)

## Groupes
K <- 3
get_group <- function(tip) {
    if (tip %in% getDescendants(tree, 105)) {
        return(2)
    }
    if (tip %in% getDescendants(tree, 110)) {
        return(3)
    }
    return(1)
}

source("./simulations/functions-anova.R")

# Computing groups
phylo_groups <- as.factor(sapply(1:n, get_group))
non_phylo_groups <- as.factor(sample(c(1, 2, 3), n, replace = TRUE))

# N répétitions pour les 2 groupes générés
phylo_groups_results <- do.call("rbind", lapply(1:N, function(id) simulate_ANOVAs(sim_id = id, groups = phylo_groups, tree = tree)))
non_phylo_groups_results <- do.call("rbind", lapply(1:N, function(id) simulate_ANOVAs(sim_id = id, groups = non_phylo_groups, tree = tree)))

# TODO : Regarder la notice de lmertest pour l'implémentation de Satterthwaite
# TODO : En utilisant l'arbre étoile, on obtient un modele mixte classique donc on peut appliquer lmerTest