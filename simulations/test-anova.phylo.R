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

# Arbre
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

# Computing groups
phylo_groups <- as.factor(sapply(1:n, get_group))
non_phylo_groups <- as.factor(sample(c(1, 2, 3), n, replace = TRUE))

overall_p <- function(my_model) {
    f <- summary(my_model)$fstatistic
    p <- pf(f[1], f[2], f[3], lower.tail = F)
    attributes(p) <- NULL
    return(p)
}

compute_y <- function(mu_vect, groups) {
    rowSums(sapply(seq_along(mu_vect), function(i) mu_vect[i] * (groups == i)))
}

# TODO : Regarder correspondance OU avec MB(+erreur de mesures)
# TODO : Refaire avec un Ornhstein-Uhlenbeck

# Code for one simulation
simulate_ <- function(sim_id,
    groups,
    n = 100,
    stoch_process = "BM",
    tree = tree,
    mu_vect = c(2, -5, 2),
    risk_threshold = 0.05,
    sub_branches = 0,
    sigma2_measure_err = 1,
    sigma2_intra_species = 1) {



    # What hypo are we testing ?
    is_H0 <- length(unique(mu_vect)) == 1

    # Are we adding sub-branches ?
    if (sub_branches != 0) {
        ## TODO: rajouter 3 petites branches au bout de l'arbre pour illustrer la variabilité intra-espece.
        ## regarder si ça dégrade la performance
        # TODO: Add sub-branching
        stop("The sub branches needs to be implemented.")
    }


    # Continuous phylo trait
    trait <- rTrait(1, tree, stoch_process)

    # Adding measure noise to the trait
    trait <- trait + rnorm(n, mean = 0, sqrt(sigma2_measure_err))

    # Simulation
    ## Réponse
    y <- compute_y(mu_vect = mu_vect, groups)
    y <- y + trait

    ## ANOVAs
    fit_ANOVA <- lm(y ~ groups)
    fitphy_ANOVA <- phylolm(y ~ groups, phy = tree, model = stoch_process)

    ## TODO refaire avec ces modalités et évaluer les erreurs de type 1 et erreurs de type 2
    ## faire scénario H_0: mu egaux -> ANOVA se plante car dep entre les indivs
    ## faire scenario H_1: mu differents -> supp ANOVA phylo se plante car pas de dep entre indiv

    methods <- as.factor(c("ANOVA", "ANOVA-Phylo"))
    
    if(is_H0){
        correct_hypothesis <- rep("H0", 2)
        
        has_selected_correctly <- c(
            overall_p(summary(fit_ANOVA)) > risk_threshold,
            overall_p(summary(fitphy_ANOVA)) > risk_threshold
        )
    } else {
        correct_hypothesis <- rep("H1", 2)

        # If the p_value is below the risk_threshold the H0 is rejected
        has_selected_correctly <- c(
            overall_p(summary(fit_ANOVA)) <= risk_threshold,
            overall_p(summary(fitphy_ANOVA)) <= risk_threshold
        )
    }

    results <- data.frame(
        sim_id = rep(sim_id, 2),
        methods = methods,
        correct_hypothesis = correct_hypothesis,
        has_selected_correctly = has_selected_correctly
    )

    return(results)
}



# TODO : Regarder la notice de lmertest pour l'implémentation de Satterthwaite
# TODO : En utilisant l'arbre étoile, on obtient un modele mixte classique donc on peut appliquer lmerTest