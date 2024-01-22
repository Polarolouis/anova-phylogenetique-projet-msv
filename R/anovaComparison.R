# Phylocomparison tools
library(phylolm)
library(phylotools)
library(phytools)
library(phylolimma)
library(ape)
library(tidyverse)

# Plotting
library(ggplot2)

# Sourcing the utils
source("./R/utils.R")

# Fixing randomness for reproducibility
set.seed(1234)

# Parameters
nb_species <- 100

# Generating the phylo tree
tree <- rphylo(nb_species, birth = 0.1, death = 0)

# Group selections tries to have group of same size
plotTree(tree, node.numbers = TRUE)

# Here I chose two ancestors to split in two the tree
ancestors <- c(102, 104)
K <- length(ancestors) # The number of groups

# I assign the groups numbers
## Matching the phylogeny
phylo_matching_groups <- sapply(1:nb_species, function(tip) {
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
    plot(tree, show.tip.label = FALSE, x.lim = 50)
    tiplabels(bg = groups, pch = 21)
    text(x = 10, y = 0, label = "This tree will be normalised.")
}

# Saving trees
png(file = "img/group_phylo_matching_tree.png")
plot_group_on_tree(tree, group = phylo_matching_groups)
dev.off()

png(file = "img/group_random_tree.png")
plot_group_on_tree(tree, group = random_groups)
dev.off()

# Normalising tree edge length
taille_tree <- diag(vcv(tree))[1]
tree$edge.length <- tree$edge.length / taille_tree

#' Returns pvalues for both F test for anova and phylogenetic anova
#'
#' @description
# TODO Describe
phyloanova_anova_pvalues <- function(
    traits, groups, tree, stoch_process,
    test_method, measurement_error = TRUE) {
    # For phylo matching
    anova_res <- lm(traits ~ groups)

    # TODO Handle the stoch process and model for phylolm (OU etc)
    model <- stoch_process

    phyloanova_res <- phylolm(traits ~ groups,
        phy = tree,
        model = model,
        measurement_error = measurement_error # To let phylolm know if there's measurement error
    )

    anova_res <- lm(traits ~ groups)
    anova_F_stat <- summary(anova_res)$fstatistic[1]
    anova_df1 <- summary(anova_res)$fstatistic[2]
    anova_df2 <- summary(anova_res)$fstatistic[3]
    anova_p_value <- pvalue_F_test(anova_F_stat,
        df1 = anova_df1, df2 = anova_df2
    )

    if (test_method %in% c("vanilla", "satterthwaite")) {
        phyloanova_F_stat <- compute_F_statistic(
            r_squared = phyloanova_res$r.squared,
            df1 = K - 1,
            df2 = nb_species - K
        )

        df1 <- K - 1
        df2 <- nb_species - K

        if (test_method == "satterthwaite") {
            # For satterthwaite ddf computation
            phyloanova_res$REML <- FALSE
            df2 <- phylolimma:::ddf_satterthwaite(phyloanova_res, tree)
        }

        phyloanova_p_value <- pvalue_F_test(phyloanova_F_stat, df1 = df1, df2 = df2)
    }

    if (test_method == "likelihood_ratio") {
        # TODO Change method name to be less deceptive
        # How to obtain the loglikehood under H0 ?
        # DONE Find the correct way to do this
        # I assume that under H0 this is like saying everyone is from the same group
        h0_phyloanova <- phylolm(traits ~ 1,
            phy = tree,
            model = model,
            measurement_error = measurement_error # To let phylolm know if there's measurement error
        )
        # But this gives a LAPACK error, the system is not inversible.

        lambda_ratio_stat <- -2*(h0_phyloanova$logLik - phyloanova_res$logLik)


        # Computes the pvalue from the statistic
        # df1 = K - 1
        phyloanova_p_value <- 1 - pchisq(lambda_ratio_stat, K-1)
    }

    list(
        phyloanova_p_value = phyloanova_p_value,
        anova_p_value = anova_p_value
    )
}

simulate_matching_and_random <- function(
    id, base_values,
    sigma2_phylo, sigma2_measure,
    stoch_process, test_method,
    risk_threshold = 0.05,
    correct_hypothesis = "H1") {
    matching_phylo_traits <- compute_trait_values(
        groups = phylo_matching_groups,
        base_values = base_values, tree,
        sigma2_phylo = sigma2_phylo, sigma2_measure = sigma2_measure,
        stoch_process = stoch_process
    )

    matching_pvalues <- phyloanova_anova_pvalues(
        traits = matching_phylo_traits,
        groups = phylo_matching_groups, tree, stoch_process = stoch_process,
        test_method = test_method, measurement_error = (sigma2_measure != 0)
    )

    random_groups_traits <- compute_trait_values(
        groups = random_groups,
        base_values = base_values, tree,
        sigma2_phylo = sigma2_phylo, sigma2_measure = sigma2_measure,
        stoch_process = stoch_process
    )

    random_groups_pvalues <- phyloanova_anova_pvalues(
        traits = random_groups_traits,
        groups = random_groups, tree, stoch_process = stoch_process,
        test_method = test_method, measurement_error = (sigma2_measure != 0)
    )

    # Concatenate pvalues
    pvalues <- c(unlist(matching_pvalues), unlist(random_groups_pvalues))

    if (correct_hypothesis == "H1") {
        correctly_selected <- pvalues < risk_threshold
    }
    if (correct_hypothesis == "H0"){
        correctly_selected <- pvalues >= risk_threshold
    }
    return(
        data.frame(
            sim_id = rep(id, 4),
            test_method = rep(c("phylo-anova", "anova"), 2),
            group_type = rep(c("matching", "random"), each = 2),
            pvalues = pvalues,
            correctly_selected = correctly_selected
        )
    )
}

# Parameters for the simulations
N <- 500
base_values <- c(1, 3) # The base trait to add
risk_threshold <- 0.05

sigma2_phylo <- 1
sigma2_measure <- 0
stoch_process <- "BM"
test_method <- "satterthwaite" # "vanilla" # "satterthwaite", "likelihood_ratio"
simulate_data <- function(N, base_values, risk_threshold, sigma2_phylo, sigma2_measure, stoch_process, test_method, correct_hypothesis = "H1") {
    simulated_data <- do.call("rbind", lapply(1:N, function(id) {
        simulate_matching_and_random(
            id = id, base_values = base_values,
            sigma2_phylo = sigma2_phylo, sigma2_measure = sigma2_measure,
            stoch_process = stoch_process,
            test_method = test_method,
            risk_threshold = risk_threshold,
            correct_hypothesis = "H1"
        )
    }))

    parameters <- paste0(
        " sigma2_measure = ", sigma2_measure,
        "; sigma2_phylo = ", sigma2_phylo,
        ";\nbase values = (", paste(c(base_values), collapse = ";"), ")",
        "; test method : ", test_method,
        "; correct hypothesis :", correct_hypothesis
    )

    return(list(data = simulated_data, parameters = parameters))
}

plot_data <- function(data, parameters, threshold = 0.95) {
    plot_data <- data %>%
        group_by(test_method, group_type) %>%
        summarize(power = mean(correctly_selected))

    p <- ggplot(plot_data, aes(x = test_method, y = power, fill = group_type)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_y_continuous(limits = c(0, 1)) +
        labs(
            title = paste0("Power vs Tested Method (", stoch_process, ") | N = ", N, ";", parameters),
            x = "Tested Method",
            y = "Power"
        ) +
        geom_hline(yintercept = threshold) +
        theme_minimal()
    p

    return(p)
}

# Vanilla
vanilla_results <- simulate_data(N, base_values, risk_threshold, sigma2_phylo,
    sigma2_measure, stoch_process,
    test_method = "vanilla"
)
vanilla_data <- vanilla_results$data
vanilla_parameters <- vanilla_results$parameters
plot_data(vanilla_data, vanilla_parameters)

vanilla_results_H0 <- simulate_data(N, base_values = c(1,1), risk_threshold, sigma2_phylo,
    sigma2_measure, stoch_process,
    test_method = "vanilla",
    correct_hypothesis = "H0"
)
vanilla_data_H0 <- vanilla_results_H0$data
vanilla_parameters_H0 <- vanilla_results_H0$parameters
plot_data(vanilla_data_H0, vanilla_parameters_H0, threshold = 0.05)

# Satterthwaite

satterthwaite_results <- simulate_data(N, base_values, risk_threshold, sigma2_phylo,
    sigma2_measure = 1, stoch_process,
    test_method = "satterthwaite"
)
satterthwaite_data <- satterthwaite_results$data
satterthwaite_parameters <- satterthwaite_results$parameters
plot_data(satterthwaite_data, satterthwaite_parameters)

satterthwaite_results_H0 <- simulate_data(N, base_values = c(1,1), risk_threshold, sigma2_phylo,
    sigma2_measure = 1, stoch_process,
    test_method = "satterthwaite", correct_hypothesis = "H0"
)
satterthwaite_data_H0 <- satterthwaite_results_H0$data
satterthwaite_parameters_H0 <- satterthwaite_results_H0$parameters
plot_data(satterthwaite_data_H0, satterthwaite_parameters_H0, threshold = 0.05)

# Likelihood ratio

lik_results <- simulate_data(N, base_values, risk_threshold, sigma2_phylo,
    sigma2_measure, stoch_process,
    test_method = "likelihood_ratio"
)
lik_data <- lik_results$data
lik_parameters <- lik_results$parameters
plot_data(lik_data, lik_parameters)

# TODO Adapt to the current code
# ## Standardized parameters
# total_variance <- 1.0 # sigma2_phylo + sigma2_error, fixed [as tree_height = 1]
# heri <- c(0.0, 0.5, 1.0) # heritability her = sigma2_phylo / total_variance. 0 means only noise. 1 means only phylo.
# snr <- 1 # signal to noise ratio snr = size_effect / total_variance

# ## Try several parameter values
# for (her in heri) {
#   res_sim <- plot_different_sigmas(sigma2_measure_err = (1 - her) * total_variance,
#                                    sigma2_intra_species = her * total_variance,
#                                    mu_vect_different = c(0, snr * total_variance, -snr * total_variance))
#   res_sim_plot <- res_sim$plot
#   res_sim_plot
#   ggsave(paste0("img/simulation_power_BM_her_", her, ".png"), plot = res_sim_plot)
# }