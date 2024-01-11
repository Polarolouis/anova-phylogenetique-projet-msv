library(phylolm)
library(phylotools)
library(phytools)
library(ape)
library(ggplot2)
library(gridExtra)
library(grid)
set.seed(1234)

N <- 500 # Number of different simulations
n <- 100

# Arbre
tree <- rphylo(n, 0.1, 0)

# Groupes
K <- 3
get_group <- function(tip) {
    if (tip %in% getDescendants(tree, 107)) {
        return(2)
    }
    if (tip %in% getDescendants(tree, 111)) {
        return(3)
    }
    return(1)
}


# Group on trees
plot_group_on_tree <- function(tree, group) {
    plot(tree, show.tip.label = FALSE, x.lim = 50)
    tiplabels(bg = group, pch = 21)
}

source("./simulations/functions-anova.R")

# Computing groups
phylo_groups <- as.factor(sapply(1:n, get_group))
non_phylo_groups <- as.factor(sample(c(1, 2, 3), n, replace = TRUE))

# Saving trees
png(file = "img/group_phylo_tree.png")
plot_group_on_tree(tree, group = phylo_groups)
dev.off()

png(file = "img/group_nonphylo_tree.png")
plot_group_on_tree(tree, group = non_phylo_groups)
dev.off()

# Normalising tree
taille_tree <- diag(vcv(tree))[1]
tree$edge.length <- tree$edge.length / taille_tree

calcul_proportion_bonne_selection <- function(data, test_method) {
    # Si mu_type = different, alors proportion de H_1 sous (H_1), donc \beta (1-\beta = puissance)
    # Si mu_type = equals, alors je calcule la proportion de H_0 bien détectés sous H_0 (vrai)
    mean(data[which(data$tested_method == test_method), ]$has_selected_correctly)
}

plot_different_sigmas <- function(sigma2_measure_err, sigma2_intra_species, mu_vect_different = c(1,4)) {
    # Mu tous différents
    mu_vect <- mu_vect_different
    # N répétitions pour les 2 groupes générés
    mu_diff_phylo_groups_results <- do.call("rbind", lapply(1:N, function(id) {
        simulate_ANOVAs(
            sim_id = id, groups = phylo_groups, tree = tree, n = n,
            df1 = K-1, df2 = n-K, 
            mu_vect = mu_vect, sigma2_measure_err = sigma2_measure_err,
            sigma2_intra_species = sigma2_intra_species
        )
    }))
    mu_diff_non_phylo_groups_results <- do.call("rbind", lapply(1:N, function(id) {
        simulate_ANOVAs(
            sim_id = id,
            groups = non_phylo_groups, 
            tree = tree, n = n, mu_vect = mu_vect,
            df1 = K-1, df2 = n-K,
            sigma2_measure_err = sigma2_measure_err, 
            sigma2_intra_species = sigma2_intra_species
        )
    }))

    puissance_mu_diff_phylo_for_phylo_groups <- calcul_proportion_bonne_selection(mu_diff_phylo_groups_results, "ANOVA-Phylo")
    puissance_mu_diff_classic_for_phylo_groups <- calcul_proportion_bonne_selection(mu_diff_phylo_groups_results, "ANOVA")

    puissance_mu_diff_phylo_for_non_phylo_groups <- calcul_proportion_bonne_selection(mu_diff_non_phylo_groups_results, "ANOVA-Phylo")
    puissance_mu_diff_classic_for_non_phylo_groups <- calcul_proportion_bonne_selection(mu_diff_non_phylo_groups_results, "ANOVA")

    # Mu égaux
    mu_vect <- rep(1, K)
    # N répétitions pour les 2 groupes générés
    mu_equals_phylo_groups_results <- do.call(
        "rbind",
        lapply(
            1:N,
            function(id) {
                simulate_ANOVAs(
                    sim_id = id, groups = phylo_groups,
                    df1 = K-1, df2 = n-K,
                    tree = tree, n = n, mu_vect = mu_vect, sigma2_measure_err = sigma2_measure_err, sigma2_intra_species = sigma2_intra_species
                )
            }
        )
    )
    mu_equals_non_phylo_groups_results <- do.call(
        "rbind",
        lapply(
            1:N, 
        function(id) {
            simulate_ANOVAs(
                sim_id = id,
                groups = non_phylo_groups, 
                tree = tree, n = n, mu_vect = mu_vect,
                df1 = K-1, df2 = n-K,
                sigma2_measure_err = sigma2_measure_err, sigma2_intra_species = sigma2_intra_species
            )
        })
    )

    # Calcul de puissance
    puissance_mu_equals_phylo_for_phylo_groups <- calcul_proportion_bonne_selection(mu_equals_phylo_groups_results, "ANOVA-Phylo")
    puissance_mu_equals_classic_for_phylo_groups <- calcul_proportion_bonne_selection(mu_equals_phylo_groups_results, "ANOVA")

    puissance_mu_equals_phylo_for_non_phylo_groups <- calcul_proportion_bonne_selection(mu_equals_non_phylo_groups_results, "ANOVA-Phylo")
    puissance_mu_equals_classic_for_non_phylo_groups <- calcul_proportion_bonne_selection(mu_equals_non_phylo_groups_results, "ANOVA")

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

    p <- ggplot(plot_data, aes(x = tested_method, y = puissance, fill = group_type)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_y_continuous(limits = c(0, 1)) +
        labs(
            title = paste0("Proportion de bonnes sélections vs Tested Method (BM) | N = ", N, "; sigma2_measure_err = ", sigma2_measure_err, "; sigma2_intra = ", sigma2_intra_species),
            x = "Tested Method",
            y = "Proportion de bonnes sélections"
        ) +
        geom_hline(yintercept = 0.95) +
            facet_grid(cols = vars(mu_type)) +
            theme_minimal()
        
    print(mu_vect_different)
    return(list(plot = p,
        mu_diff_phylo_groups_results = mu_diff_phylo_groups_results,
        mu_equals_phylo_groups_results = mu_equals_phylo_groups_results,
        mu_diff_non_phylo_groups_results = mu_diff_non_phylo_groups_results,
        mu_equals_non_phylo_groups_results = mu_equals_non_phylo_groups_results
    ))
}
# TODO : Regarder la notice de lmertest pour l'implémentation de Satterthwaite
# TODO : En utilisant l'arbre étoile, on obtient un modele mixte classique donc on peut appliquer lmerTest
ggsave <- function(..., bg = "white") ggplot2::ggsave(..., bg = bg)

# Seulement l'information phylo
only_phylo_data <- plot_different_sigmas(0, 1, mu_vect_different = c(1,5,10))
plot_only_phylo <- only_phylo_data$plot
plot_only_phylo
# 2 erreurs
both_errors_identical_data <- plot_different_sigmas(1, 1, mu_vect_different = c(1, 5, 10))
plot_both_errors_identical <- both_errors_identical_data$plot
plot_both_errors_identical

ggsave("img/simulation_power_pure_BM.png", plot = plot_only_phylo) # BM pur
ggsave("img/simulation_power_pure_both_error.png", plot = plot_both_errors_identical) # Both error same
# ggsave("img/simulation_power_pure_both_error_more_phylo.png", plot = plot_different_sigmas(1, 2)) # Both error same
