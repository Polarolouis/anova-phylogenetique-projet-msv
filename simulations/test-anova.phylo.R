library(phylolm)
library(phylotools)
library(phytools)
library(ape)
set.seed(1234)

N <- 100 # Number of different simulations

# Preparing output lists
# simulation_results <- data.frame(
#     sim_id = numeric(),
#     positive_classic_r_squared = numeric(),
#     positive_phylo_r_squared = numeric(),
#     positive_classic_adjusted_r_squared = numeric(),
#     positive_phylo_adjusted_r_squared = numeric(),
#     negative_classic_r_squared = numeric(),
#     negative_phylo_r_squared = numeric(),
#     negative_classic_adjusted_r_squared = numeric(),
#     negative_phylo_adjusted_r_squared = numeric(),
#     row.names = NULL, check.rows = FALSE, check.names = TRUE,
#     stringsAsFactors = default.stringsAsFactors()
# )

tree <- rphylo(n, 0.1, 0)

## Groupes
get_group <- function(tip) {
    if (tip %in% getDescendants(tree, 105)) {
        return(2)
    }
    if (tip %in% getDescendants(tree, 110)) {
        return(3)
    }
    return(1)
}

phylo_group <- as.factor(sapply(1:n, get_group))

# Code for one simulation
simulate_positive_negative <- function(sim_id, n = 100, stoch_process = "BM", tree = tree, phylo_group = phylo_group) {

    sigma2err <- 1

    # Continuous phylo trait
    trait <- rTrait(1, tree, stoch_process)

    # Adding noise to the trait
    trait <- trait + rnorm(n, mean = 0, sqrt(sigma2err))

    # Simulation positive
    # TODO : Refaire avec un Ornhstein-Uhlenbeck
        ## TODO refaire avec ces modalités et évaluer les erreurs de type 1 et erreurs de type 2
        ## faire scénario H_0: mu egaux -> ANOVA se plante car dep entre les indivs
        ## faire scenario H_1: mu differents -> supp ANOVA phylo se plante car pas de dep entre indiv
        ## TODO: rajouter 3 petites branches au bout de l'arbre pour illustrer la variabilité intra-espece.
        ## regarder si ça dégrade la performance
    # TODO : Regarder correspondance OU avec MB(+erreur de mesures)

    ## Réponse
    mu1 <- 2
    mu2 <- -5
    mu3 <- 2
    y <- mu1 * (phylo_group == 1) + mu2 * (phylo_group == 2) + mu3 * (phylo_group == 3)
    y <- y + trait

    # par(mar = c(5, 0, 0, 0) + 0.1)
    # plot(tree, show.tip.label = FALSE, x.lim = 50)
    # tiplabels(bg = group, pch = 21)
    # phydataplot(y, tree, scaling = 0.1, offset = 4)

    pos_fit_ANOVA <- lm(y ~ phylo_group)

    pos_fitphy_ANOVA <- phylolm(y ~ phylo_group, phy = tree)

    # Simulation négative

    groups_non_phylo <- as.factor(sample(c(1, 2, 3), n, replace = TRUE))
    y_non_phy <- mu1 * (groups_non_phylo == 1) + mu2 * (groups_non_phylo == 2) + mu3 * (groups_non_phylo == 3)
    y_non_phy <- y_non_phy + trait

    # par(mar = c(5, 0, 0, 0) + 0.1)
    # plot(tree, show.tip.label = FALSE, x.lim = 50)
    # tiplabels(bg = groups_non_phylo, pch = 21)
    # phydataplot(y_non_phy, tree, scaling = 0.1, offset = 4)

    neg_fit_ANOVA <- lm(y_non_phy ~ groups_non_phylo)
    neg_fitphy_ANOVA <- phylolm(y_non_phy ~ groups_non_phylo, phy = tree)

    # Summary
    ## Positive
    pos_fit_summary <- summary(pos_fit_ANOVA)
    pos_fitphy_summary <- summary(pos_fitphy_ANOVA)

    ## Negative
    neg_fit_summary <- summary(neg_fit_ANOVA)
    neg_fitphy_summary <- summary(neg_fitphy_ANOVA)
    return(data.frame(
        sim_id = sim_id,
        positive_classic_r_squared = pos_fit_summary$r.squared,
        positive_phylo_r_squared = pos_fitphy_summary$r.squared,
        positive_classic_adjusted_r_squared = pos_fit_summary$adj.r.squared,
        positive_phylo_adjusted_r_squared = pos_fitphy_summary$adj.r.squared,
        negative_classic_r_squared = neg_fit_summary$r.squared,
        negative_phylo_r_squared = neg_fitphy_summary$r.squared,
        negative_classic_adjusted_r_squared = neg_fit_summary$adj.r.squared,
        negative_phylo_adjusted_r_squared = neg_fitphy_summary$adj.r.squared,
        row.names = NULL, check.rows = FALSE, check.names = TRUE,
        stringsAsFactors = default.stringsAsFactors()
    ))
}


simulation_results <- do.call(rbind, lapply(seq(N), function(sim_id) {
    simulate_positive_negative(sim_id)
}))

# TODO : Regarder la notice de lmertest pour l'implémentation de Satterthwaite
# TODO : En utilisant l'arbre étoile, on obtient un modele mixte classique donc on peut appliquer lmerTest