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
nb_species <- 20

# Generating the phylo tree
tree <- rphylo(nb_species, birth = 0.1, death = 0)

# Group selections tries to have group of same size
plotTree(tree, node.numbers = TRUE)

# Here I chose two ancestors to split in two the tree
ancestors <- c(22, 23)
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


# TODO REMOVE
# simulate_matching_and_random <- function(
#     id, base_values,
#     sigma2_phylo, sigma2_measure,
#     stoch_process, test_method,
#     risk_threshold = 0.05,
#     correct_hypothesis = "H1") {
#     # To be in the right configuration for typeI error
#     if (correct_hypothesis == "H0") {
#         base_values <- rep(0, length(base_values))
#     }
#     matching_phylo_traits <- compute_trait_values(
#         groups = phylo_matching_groups,
#         base_values = base_values, tree,
#         sigma2_phylo = sigma2_phylo, sigma2_measure = sigma2_measure,
#         stoch_process = stoch_process
#     )
#     matching_pval_df <- phyloanova_anova_pvalues(
#         traits = matching_phylo_traits,
#         groups = phylo_matching_groups, tree, stoch_process = stoch_process,
#         test_method = test_method, measurement_error = TRUE
#     )
#     matching_pvalues <- matching_pval_df[c(1, 2)]

#     matching_df2 <- matching_pval_df[c(3, 4)]

#     random_groups_traits <- compute_trait_values(
#         groups = random_groups,
#         base_values = base_values, tree,
#         sigma2_phylo = sigma2_phylo, sigma2_measure = sigma2_measure,
#         stoch_process = stoch_process
#     )

#     random_groups_pval_df2 <- phyloanova_anova_pvalues(
#         traits = random_groups_traits,
#         groups = random_groups, tree, stoch_process = stoch_process,
#         test_method = test_method, measurement_error = TRUE
#     )
#     random_groups_pvalues <- random_groups_pval_df2[c(1, 2)]

#     random_groups_df2 <- random_groups_pval_df2[c(3, 4)]
#     # Concatenate pvalues
#     pvalues <- c(unlist(matching_pvalues), unlist(random_groups_pvalues))

#     if (correct_hypothesis == "H1") {
#         correctly_selected <- pvalues < risk_threshold
#     }
#     if (correct_hypothesis == "H0") {
#         correctly_selected <- pvalues >= risk_threshold
#     }
#     return(
#         data.frame(
#             sim_id = rep(id, 4),
#             anova_model = rep(c("phylo-anova", "anova"), 2),
#             group_type = rep(c("matching", "random"), each = 2),
#             pvalues = pvalues,
#             correctly_selected = correctly_selected,
#             df2 = unlist(c(matching_df2, random_groups_df2))
#         )
#     )
# }



# compare_methods <- function(
#     N, base_values, risk_threshold, sigma2_phylo,
#     sigma2_measure, stoch_process, methods_to_test = c("vanilla", "satterthwaite"), correct_hypothesis = "H1") {
#     if (any(!(methods_to_test %in% c("vanilla", "satterthwaite", "lrt")))) {
#         stop("Unknown method to test.")
#     }
#     #  Generating data for each method
#     # TODO Utiliser les mêmes données pour les méthodes
#     ##  To compute power
#     full_power_data <-
#         do.call("rbind", lapply(methods_to_test, function(method) {
#             data <- simulate_data(
#                 N = N,
#                 base_values = base_values,
#                 risk_threshold = risk_threshold,
#                 sigma2_phylo = sigma2_phylo,
#                 sigma2_measure = sigma2_measure,
#                 test_method = method,
#                 stoch_process = stoch_process,
#                 correct_hypothesis = "H1"
#             )$data
#             #  Adding a column to identify the approximation method
#             data$tested_method <- method
#             data$metric_type <- "power"
#             data
#         }))
#     ##  To compute type I error
#     full_typeI_data <-
#         do.call("rbind", lapply(methods_to_test, function(method) {
#             data <- simulate_data(
#                 N = N,
#                 base_values = base_values,
#                 risk_threshold = risk_threshold,
#                 sigma2_phylo = sigma2_phylo,
#                 sigma2_measure = sigma2_measure,
#                 test_method = method,
#                 stoch_process = stoch_process,
#                 correct_hypothesis = "H0"
#             )$data
#             #  Adding a column to identify the approximation method
#             data$tested_method <- method
#             data$metric_type <- "typeI"
#             data
#         }))

#     data <- rbind(full_power_data, full_typeI_data)
#     return(
#         list(
#             data = data,
#             sim_parameters = list(
#                 N = N,
#                 base_values = base_values,
#                 risk_threshold = risk_threshold,
#                 sigma2_phylo = sigma2_phylo,
#                 sigma2_measure = sigma2_measure,
#                 stoch_process = stoch_process
#             )
#         )
#     )
# }

# plot_simulation_data <- function(data, parameters_string, threshold = 0.95) {
#     plot_data <- data %>%
#         group_by(anova_model, group_type) %>%
#         summarize(power = mean(correctly_selected))

#     p <- ggplot(plot_data, aes(x = anova_model, y = power, fill = group_type)) +
#         geom_bar(stat = "identity", position = "dodge") +
#         scale_y_continuous(limits = c(0, 1)) +
#         labs(
#             title = paste0("Power vs Tested Method (", stoch_process, ") | N = ", N, ";", parameters_string),
#             x = "Tested Method",
#             y = "Power"
#         ) +
#         geom_hline(yintercept = threshold) +
#         theme_minimal()
#     p

#     return(p)
# }

plot_comparison <- function(data, sim_parameters) {
    #  Retrieving simulation parameters
    risk_threshold <- sim_parameters$risk_threshold
    N <- sim_parameters$N
    sigma2_measure <- sim_parameters$sigma2_measure
    sigma2_phylo <- sim_parameters$sigma2_phylo
    base_values <- sim_parameters$base_values
    stoch_process <- sim_parameters$stoch_process

    plot_title <- paste0(
        "N = ", N, ";", " sigma2_measure = ", sigma2_measure,
        "; sigma2_phylo = ", sigma2_phylo,
        ";\nbase values = (", paste(c(base_values), collapse = ","), ");",
        "\nStoch process : ", stoch_process
    )

    #  Preparing plot data
    plot_data <- data %>%
        group_by(tested_method, anova_model, group_type, metric_type) %>%
        summarize(metric = mean(correctly_selected))
    #  Reversing the metric to really be typeI error (ie the prop of errors made)
    plot_data[plot_data$metric_type == "typeI", ] <- plot_data[plot_data$metric_type == "typeI", ] %>% mutate(metric = 1 - metric)

    # Adding a threshold
    plot_data <- plot_data %>%
        ungroup() %>%
        mutate(
            hline_risk_threshold = case_when(
                plot_data$metric_type == "power" ~ -0.1,
                plot_data$metric_type == "typeI" ~ risk_threshold
            )
        )
    #  To be out of bounds

    p <- ggplot(plot_data, aes(x = anova_model, y = metric, fill = group_type)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_text(aes(label = round(metric, digits = 3)), vjust = -0.5, position = position_dodge(width = 0.9)) +
        scale_y_continuous(limits = c(0, 1.2)) +
        labs(
            title = plot_title,
            x = "Anova Method",
            y = "Metric"
        ) +
        theme_minimal()
    p <- p + facet_grid(tested_method ~ metric_type)
    p <- p + geom_hline(aes(yintercept = hline_risk_threshold))

    return(p)
}

#  Comparing methods
# comparison <- compare_methods(N,
#     base_values = c(1, 2), risk_threshold, sigma2_phylo = 0,
#     sigma2_measure = 0.5, stoch_process, methods_to_test = c("vanilla", "satterthwaite", "lrt")
# )
# comparison_data <- comparison$data

# plot_comparison(comparison_data, sim_parameters = comparison$sim_parameters)


## Standardized parameters
total_variance <- 1.0 # sigma2_phylo + sigma2_error, fixed [as tree_height = 1]
heri <- c(0.0, 0.25, 0.5, 1.0) # heritability her = sigma2_phylo / total_variance. 0 means only noise. 1 means only phylo.
snr <- 1 # signal to noise ratio snr = size_effect / total_variance

## Try several parameter values
ggsave <- function(..., bg = "white") ggplot2::ggsave(..., bg = bg)
for (her in heri) {
    sim <- compare_methods(N,
        base_values = c(0, snr * total_variance), risk_threshold, sigma2_phylo = her * total_variance,
        sigma2_measure = (1 - her) * total_variance, stoch_process, methods_to_test = c("vanilla", "satterthwaite", "lrt")
    )

    res_sim_plot <- plot_comparison(sim$data, sim$sim_parameters)
    res_sim_plot
    ggsave(paste0("img/simulation_BM_her_", her, ".png"), plot = res_sim_plot)
}
