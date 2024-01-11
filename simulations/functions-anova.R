overall_p <- function(my_model) {
    f <- summary(my_model)$fstatistic
    p <- pf(f[1], f[2], f[3], lower.tail = F)
    attributes(p) <- NULL
    return(p)
}

compute_F_statistic <- function(r_squared, df1, df2) {
    # df1 = k, le nombre de prédicteur
    # df2 = n - (k+1), n le nombre d'observation
    return(r_squared / (1 - r_squared) * df2 / df1)
}

phylo_p_value <- function(r_squared, df1, df2) {
    F_stat <- compute_F_statistic(r_squared, df1, df2)
    return(1 - pf(F_stat, df1 = df1, df2 = df2))
}

compute_y <- function(mu_vect, groups) {
    rowSums(sapply(seq_along(mu_vect), function(i) mu_vect[i] * (groups == i)))
}

# TODO : Regarder correspondance OU avec MB(+erreur de mesures)
# TODO : Refaire avec un Ornhstein-Uhlenbeck

# Code for one simulation
simulate_ANOVAs <- function(
    sim_id,
    groups,
    tree,
    n,
    mu_vect,
    stoch_process = "BM",
    risk_threshold = 0.05,
    sub_branches = 0,
    sigma2_measure_err = 1,
    sigma2_intra_species = 1,
    df1 = K - 1,
    df2 = n - K) {
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
    trait <- rTrait(1, tree, stoch_process, parameters = list(sigma2 = sigma2_intra_species))

    # Adding measure noise to the trait
    trait <- trait + rnorm(n, mean = 0, sqrt(sigma2_measure_err))

    # Simulation
    ## Réponse
    y <- compute_y(mu_vect = mu_vect, groups)
    y <- y + trait

    ## ANOVAs
    fit_ANOVA <- lm(y ~ groups)
    if (stoch_process == "OU") {
        model = "OUfixedRoot"
    } else {
        model = stoch_process
    }
    fitphy_ANOVA <- phylolm(y ~ groups, phy = tree, model = model, measurement_error = (sigma2_measure_err != 0))

    ## DONE refaire avec ces modalités et évaluer les erreurs de type 1 et erreurs de type 2
    ## faire scénario H_0: mu egaux -> ANOVA se plante car dep entre les indivs
    ## faire scenario H_1: mu differents -> supp ANOVA phylo se plante car pas de dep entre indiv

    tested_method <- as.factor(c("ANOVA", "ANOVA-Phylo"))

    if (is_H0) {
        correct_hypothesis <- rep("H0", 2)

        has_selected_correctly <- c(
            overall_p(fit_ANOVA) > risk_threshold,
            phylo_p_value(fitphy_ANOVA$r.squared, df1, df2) > risk_threshold
        )
        selected_hypothesis <- sapply(1:2, function(id) {
            if (has_selected_correctly[id]) {
                return("H0")
            } else {
                return("H1")
            }
        })
    } else {
        correct_hypothesis <- rep("H1", 2)

        # If the p_value is below the risk_threshold the H0 is rejected
        has_selected_correctly <- c(
            overall_p(fit_ANOVA) <= risk_threshold,
            phylo_p_value(fitphy_ANOVA$r.squared, df1, df2) <= risk_threshold
        )
        selected_hypothesis <- sapply(1:2, function(id) {
            if (has_selected_correctly[id]) {
                return("H1")
            } else {
                return("H0")
            }
        })
    }

    results <- data.frame(
        sim_id = rep(sim_id, 2),
        tested_method = tested_method,
        correct_hypothesis = correct_hypothesis,
        selected_hypothesis = selected_hypothesis,
        has_selected_correctly = has_selected_correctly
    )

    return(results)
}