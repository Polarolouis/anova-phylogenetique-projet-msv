#' Get the group number of a tip of the given tree
#'
#' @description
#' Returns the group number (based on the number of ancestors provided)
#' for the given tip
get_phylo_group <- function(tip, tree, ancestors = c(102, 104)) {
    # Sanity checks
    if (!all(ancestors %in% tree$edge[, 1])) {
        stop("At least one ancestor is not in the given tree")
    }
    if (!paste0("t", tip) %in% tree$tip.label) {
        stop("Provided tip is not in the tree.")
    }

    for (ancestor_id in seq_along(ancestors)) {
        if (tip %in% getDescendants(tree, ancestors[ancestor_id])) {
            # If the tip is a descendant of this ancestor return its index
            return(ancestor_id)
        }
    }

    # If we reach this line the tip is not part of the given ancestors
    warning(paste0(
        "The tip ", tip,
        " is not a descendant of the provided ancestors"
    ))
}

#' Computes the F_statistic from the r_squared
#'
#' @description
#' Use the formula between r_squared and df1 (K-1), df2 (n - K) to return the
#' F statistic to test against the Fisher distribution.
compute_F_statistic <- function(r_squared, df1, df2) {
    # df1 = k - 1, le nombre de prédicteur
    # df2 = n - k, n le nombre d'observation
    return(r_squared / (1 - r_squared) * df2 / df1)
}

#' Computes the pvalue from an F statistic
pvalue_F_test <- function(F_stat, df1, df2) {
    return(unname(1 - pf(F_stat, df1, df2)))
}

#' @title Get the number of species
#'
#' @description
#' Compute the number of different species on a tree that possibly has replicates
#' coded as tips with zero length branches.
#'
#' From phylolimma
#'
#' @param phy a phylogentic tree, with possible replicates coded as tips with zero length branches.
#' @param tol a numeric value giving the tolerance to consider a branch length significantly greater than zero.
#'
#' @return the number of different species in the tree
#'
#' @keywords internal
#'
getSpeciesNumber <- function(phy, tol = .Machine$double.eps^(1 / 2)) {
    n <- length(phy$tip.label)
    R <- ape::cophenetic.phylo(phy) <= tol
    R <- colSums(R)
    nspecies <- 0
    ind <- 1
    while (ind <= length(R)) {
        nspecies <- nspecies + 1
        ind <- ind + R[ind]
    }
    return(nspecies)
}

#' @title Function for species ddf
#'
#' From phylolimma (pbastide/phylolimma)
#'
#' @param fitlm a phylolm fit
#' @param phylo the corresponding phylogenetic tree
#'
#' @return nspecies - nvariables
#'
#' @keywords internal
#'
ddf_species <- function(fitlm, phylo) {
    nspecies <- getSpeciesNumber(phylo)
    return(list(ddf = nspecies - fitlm$d))
}

#' This code computes the satterthwaite approximation to obtain degrees of
#' freedom for a given tree
ddf_satterthwaite_sum <- function(fit_phylolm, phylo, REML = FALSE) {
    if (is.null(fit_phylolm$sigma2_error) || fit_phylolm$sigma2_error == 0) {
        #  In this case we return the classical df
        return(ddf_species(fit_phylolm, phylo))
    }

    X <- fit_phylolm$X
    y <- as.matrix(fit_phylolm$y)
    rownames(y) <- names(fit_phylolm$y)
    yhat <- as.matrix(fit_phylolm$fitted.values)
    rownames(yhat) <- names(fit_phylolm$fitted.values)
    n <- length(phylo$tip.label)
    d <- ncol(X)

    ## Likelihood function
    minusLogLik <- function(pars, y, yhat, X, phy, model) {
        n <- nrow(X)
        d <- ncol(X)
        parameters <- list(sigma2 = exp(pars[1]), sigma2_error = exp(pars[2] - pars[1]))
        phytrans <- transf.branch.lengths(phy, model, parameters = parameters)$tree
        comp <- three.point.compute(phytrans, P = y - yhat, Q = X)

        if (!REML) {
            n2llh <- as.numeric(n * log(2 * pi) + n * log(parameters$sigma2) + comp$logd + comp$PP / parameters$sigma2) # -2 log-likelihood
        } else {
            # log|X'V^{-1}X|
            ldXX <- determinant(comp$QQ, logarithm = TRUE)$modulus
            n2llh <- as.numeric((n - d) * log(2 * pi) + (n - d) * log(parameters$sigma2) + comp$logd + comp$PP / parameters$sigma2 + ldXX) # -2 log-likelihood
        }

        return(n2llh / 2)
    }

    # Using the log scale so that parameters are on the entire real line
    optpars <- c(log(fit_phylolm$sigma2), log(fit_phylolm$sigma2_error))

    # all.equal(minusLogLik(optpars, y, yhat, X, phylo, "BM"),
    #           -fit_phylolm$logLik)

    ## Hessian: numerical computation
    fun <- function(x) {
        return(minusLogLik(x, y, yhat, X, phylo, "BM"))
    }
    J <- diag(c(1 / fit_phylolm$sigma2, 1 / fit_phylolm$sigma2_error))
    A <- compute_hessian(optpars = optpars, fun = fun, grad_trans = J)

    ## Gradient of f
    K <- vcv(phylo)
    Kd <- diag(diag(K))
    V <- fit_phylolm$sigma2 * K + fit_phylolm$sigma2_error * Kd
    Vinv <- chol2inv(chol(V))
    ell <- c(0, 1)

    C <- fit_phylolm$vcov
    if (!REML) C <- C * (n - d) / n
    # Cbis <- solve(t(X) %*% Vinv %*% X)
    # all.equal(C, Cbis)

    facmat <- C %*% t(X) %*% Vinv
    derfsigma2 <- t(ell) %*% facmat %*% K %*% t(facmat) %*% ell
    derfsigma2_error <- t(ell) %*% facmat %*% Kd %*% t(facmat) %*% ell
    derf <- c(derfsigma2, derfsigma2_error)

    ## Variance estimation
    varestim <- t(derf) %*% A %*% derf

    ## Satterthwaite
    ddf <- 2 * (t(ell) %*% C %*% ell)^2 / varestim

    return(list(ddf = ddf, vcov = A))
}

# Adapted from lmerTest
# https://github.com/runehaubo/lmerTestR/blob/35dc5885205d709cdc395b369b08ca2b7273cb78/R/lmer.R#L173
compute_hessian <- function(optpars, fun, grad_trans, tol = 1e-8, ...) {
    # Compute Hessian:
    h <- numDeriv::hessian(func = fun, x = optpars, ...)
    # back transformation of parameters
    h <- t(grad_trans) %*% h %*% grad_trans
    # Eigen decompose the Hessian:
    eig_h <- eigen(h, symmetric = TRUE)
    evals <- eig_h$values
    neg <- evals < -tol
    pos <- evals > tol
    zero <- evals > -tol & evals < tol
    if (sum(neg) > 0) { # negative eigenvalues
        eval_chr <- if (sum(neg) > 1) "eigenvalues" else "eigenvalue"
        evals_num <- paste(sprintf("%1.1e", evals[neg]), collapse = " ")
        warning(sprintf(
            "Model failed to converge with %d negative %s: %s",
            sum(neg), eval_chr, evals_num
        ), call. = FALSE)
    }
    # Note: we warn about negative AND zero eigenvalues:
    if (sum(zero) > 0) { # some eigenvalues are zero
        eval_chr <- if (sum(zero) > 1) "eigenvalues" else "eigenvalue"
        evals_num <- paste(sprintf("%1.1e", evals[zero]), collapse = " ")
        warning(sprintf(
            "Model may not have converged with %d %s close to zero: %s",
            sum(zero), eval_chr, evals_num
        ))
    }
    # Compute vcov(varpar):
    pos <- eig_h$values > tol
    q <- sum(pos)
    # Using the Moore-Penrose generalized inverse for h:
    h_inv <- with(eig_h, {
        vectors[, pos, drop = FALSE] %*% diag(1 / values[pos], nrow = q) %*%
            t(vectors[, pos, drop = FALSE])
    })
    return(h_inv)
}

#' Compute trait values for the given groups
#'
#' @description
#' For each value of group, the base value is matched in the order of
#' the vector (1st value for 1st group and etc).
#' Then the phylogenetic trait corresponding from the tree is computed.
#' And finally the error (gaussian centered) is computed.
#' The sum is returned
compute_trait_values <- function(
    groups,
    base_values,
    tree,
    sigma2_phylo,
    sigma2_measure,
    stoch_process = "BM") {
    # Here we compute the base values for each of the species
    trait <- rowSums(sapply(seq_along(base_values), function(i) base_values[i] * (groups == i)))

    # The phylogenetic induced value
    trait_phylo <- rTrait(1, tree,
        stoch_process,
        parameters = list(sigma2 = sigma2_phylo)
    )

    # The measure error
    trait_error <- rnorm(length(groups), mean = 0, sd = sqrt(sigma2_measure))

    return(trait + trait_phylo + trait_error)
}

#' Infere an ANOVA and a phyloanova
#'
#' @param y the vector of traits for which to fit the models
#' @param groups the groups to which fit the models
#' @param tree the phylo tree to use
#' @param stoch_process the stochastic process to use for the phylolm
#'
infere_anova_phyloanova <- function(y, groups, tree, stoch_process = "BM") {
    #  The fits
    fit_anova <- lm(y ~ groups)
    fit_phylolm <- phylolm(y ~ groups, phy = tree, model = stoch_process, measurement_error = TRUE)
    return(list(anova = fit_anova, phyloanova = fit_phylolm))
}

# TODO renvoyer pour toutes les méthodes testées
#' Return pvalues for the anova and the phyloanova
#'
#' @param fit_anova the lm fit
#' @param fit_phylolm the phylolm fit
#' @param tree the phylo tree
#' @param tested_method the method to test should be one of : "vanilla",
#' "satterthwaite", "lrt"
pvalues_from_fits <- function(
    fit_anova,
    fit_phylolm,
    tree,
    tested_method = c("vanilla", "satterthwaite", "lrt"),
    REML = FALSE) {
    #  For sanity test
    match.arg(tested_method)

    # Extracting values
    nb_species <- nrow(model.frame(fit_anova))
    K <- length(unique(model.frame(fit_anova)$groups))

    is_invalid_value <- function(value) {
        return(is.nan(value) ||
            is.null(value) ||
            is.infinite(value) ||
            value == 0)
    }

    #  The default value for the df2
    df1 <- K - 1
    df2 <- nb_species - K

    anova_F_stat <- summary(fit_anova)$fstatistic[1]
    anova_df1 <- summary(fit_anova)$fstatistic[2]
    anova_df2 <- summary(fit_anova)$fstatistic[3]
    pvalue_anova <- pvalue_F_test(anova_F_stat,
        df1 = anova_df1, df2 = anova_df2
    )

    switch(tested_method,
        "vanilla" = {
            F_stat <- compute_F_statistic(
                r_squared = fit_phylolm$r.squared,
                df1 = df1,
                df2 = df2
            )
            pvalue_phylolm <- 1 - pf(F_stat, df1, df2)
        },
        "satterthwaite" = {
            df2 <- ddf_satterthwaite_sum(fit_phylolm = fit_phylolm, phylo = tree, REML = REML)$ddf
            F_stat <- compute_F_statistic(
                r_squared = fit_phylolm$r.squared,
                df1 = df1,
                df2 = df2
            )
            pvalue_phylolm <- 1 - pf(F_stat, df1, df2)
        },
        "lrt" = {
            h0_phylolm <- phylolm(fit_phylolm$y ~ 1,
                phy = tree,
                model = fit_phylolm$model,
                measurement_error = !is_invalid_value(fit_phylolm$sigma2_error) # To let phylolm know if there's measurement error
            )
            lambda_ratio_stat <- -2 * (h0_phylolm$logLik - fit_phylolm$logLik)


            # Computes the pvalue from the statistic
            # df1 = K - 1
            pvalue_phylolm <- 1 - pchisq(lambda_ratio_stat, df1)
        }
    )
    return(data.frame(
        tested_method = tested_method,
        pvalue_anova = pvalue_anova,
        pvalue_phylolm = pvalue_phylolm,
        anova_df1 = anova_df1,
        anova_df2 = anova_df2,
        phylolm_df1 = df1,
        phylolm_df2 = df2
    ))
}

#' Indicated if the test has selected correctly
#'
#' @param correct_hypothesis the real hypothesis in "H0" or "H1"
#' @param pvalue the pvalue to test for
#' @param risk_threshold the type I error risk threshold (default: 0.05)
test_selected_correctly <- function(
    correct_hypothesis = c("H0","H1"),
    pvalue,
    risk_threshold = 0.05) {
    # Sanity check
    match.arg(correct_hypothesis)
    
    switch(correct_hypothesis,
        "H0" = {
            has_selected_correctly <- (pvalue >= risk_threshold)
        },
        "H1" = {
            has_selected_correctly <- (pvalue < risk_threshold)
        }
    )
    return(has_selected_correctly)
}