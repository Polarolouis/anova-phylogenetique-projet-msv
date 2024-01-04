library(phylolm)
library(phylotools)
library(phytools)
library(ape)
library(ggplot2)
set.seed(1234)

N <- 500 # Number of different simulations
n <- 100

# Arbre
tree <- rphylo(n, 0.1, 0)
taille_tree <- diag(vcv(tree))[1]
tree$edge.length <- tree$edge.length/taille_tree

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

calcul_puissance <- function(data, test_method) {
    mean(data[which(data$tested_method == test_method), ]$has_selected_correctly)
}

sigma2_measure_err = 1
sigma2_intra_species = 2

# Mu tous différents
mu_vect <- c(1, 5, 10)
# N répétitions pour les 2 groupes générés
mu_diff_phylo_groups_results <- do.call("rbind", lapply(1:N, function(id) simulate_ANOVAs(sim_id = id, groups = phylo_groups, tree = tree, mu_vect = mu_vect)))
mu_diff_non_phylo_groups_results <- do.call("rbind", lapply(1:N, function(id) simulate_ANOVAs(sim_id = id, groups = non_phylo_groups, tree = tree, mu_vect = mu_vect)))

puissance_mu_diff_phylo_for_phylo_groups <- calcul_puissance(mu_diff_phylo_groups_results, "ANOVA-Phylo")
puissance_mu_diff_classic_for_phylo_groups <- calcul_puissance(mu_diff_phylo_groups_results, "ANOVA")

puissance_mu_diff_phylo_for_non_phylo_groups <- calcul_puissance(mu_diff_non_phylo_groups_results, "ANOVA-Phylo")
puissance_mu_diff_classic_for_non_phylo_groups <- calcul_puissance(mu_diff_non_phylo_groups_results, "ANOVA")

# Mu égaux
mu_vect <- rep(1, K)
# N répétitions pour les 2 groupes générés
mu_equals_phylo_groups_results <- do.call("rbind", lapply(1:N, function(id) simulate_ANOVAs(sim_id = id, groups = phylo_groups, tree = tree, mu_vect = mu_vect, sigma2_measure_err = sigma2_measure_err)))
mu_equals_non_phylo_groups_results <- do.call("rbind", lapply(1:N, function(id) simulate_ANOVAs(sim_id = id, groups = non_phylo_groups, tree = tree, mu_vect = mu_vect, sigma2_measure_err = sigma2_measure_err)))

# Calcul de puissance
puissance_mu_equals_phylo_for_phylo_groups <- calcul_puissance(mu_equals_phylo_groups_results, "ANOVA-Phylo")
puissance_mu_equals_classic_for_phylo_groups <- calcul_puissance(mu_equals_phylo_groups_results, "ANOVA")

puissance_mu_equals_phylo_for_non_phylo_groups <- calcul_puissance(mu_equals_non_phylo_groups_results, "ANOVA-Phylo")
puissance_mu_equals_classic_for_non_phylo_groups <- calcul_puissance(mu_equals_non_phylo_groups_results, "ANOVA")

# Graphiques
puissances <- c(
    puissance_mu_diff_phylo_for_phylo_groups,
    puissance_mu_diff_classic_for_phylo_groups,
    puissance_mu_equals_phylo_for_phylo_groups,
    puissance_mu_equals_classic_for_phylo_groups,
    puissance_mu_diff_phylo_for_non_phylo_groups,
    puissance_mu_diff_classic_for_non_phylo_groups,
    puissance_mu_equals_phylo_for_non_phylo_groups,
    puissance_mu_equals_classic_for_non_phylo_groups
)
plot_data <- data.frame(
    puissance = puissances,
    tested_method = rep(c("ANOVA-Phylo", "ANOVA"), 4),
    group_type = rep(c("phylo", "non_phylo"), each = 4),
    mu_type = rep(rep(c("different", "equals"), each = 2), 2)
)

ggplot(plot_data, aes(x = tested_method, y = puissance, fill = interaction(group_type, mu_type))) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
        title = paste0("Proportions correctes vs Tested Method (BM) | N = ", N, "; sigma_measure = ", sigma2_measure_err, "; sigma2_intra = ", sigma2_intra_species),
        x = "Tested Method",
        y = "Proportions correctes"
    ) +
    geom_hline(yintercept = 0.95) +
    theme_minimal()
# TODO : Regarder la notice de lmertest pour l'implémentation de Satterthwaite
# TODO : En utilisant l'arbre étoile, on obtient un modele mixte classique donc on peut appliquer lmerTest

