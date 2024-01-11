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