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
phylomatching_groups <- sapply(1:nb_species, function(tip) {
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
plot_group_on_tree(tree, group = phylomatching_groups)
dev.off()

png(file = "img/group_random_tree.png")
plot_group_on_tree(tree, group = random_groups)
dev.off()

# Normalising tree edge length
taille_tree <- diag(vcv(tree))[1]
tree$edge.length <- tree$edge.length / taille_tree

#' Simulates 2 sets of data, one matching phylo, one random and apply the
#' 3 methods on it (vanilla, satterthwaite, lrt)
simulate_all_methods <- function(
    sim_id, correct_hypothesis = c("H0", "H1"), base_values,
    sigma2_phylo,
    sigma2_measure, risk_threshold = 0.05) {
    # Be sure to ignore base_values if testing H0
    if (correct_hypothesis == "H0") {
        base_values <- rep(0, length(base_values))
    }

    # Computing traits
    phylomatching_y <- compute_trait_values(
        groups = phylomatching_groups,
        base_values = base_values,
        tree = tree, sigma2_phylo = sigma2_phylo,
        sigma2_measure = sigma2_measure
    )
    random_y <- compute_trait_values(
        groups = random_groups,
        base_values = base_values,
        tree = tree, sigma2_phylo = sigma2_phylo,
        sigma2_measure = sigma2_measure
    )

    # Fits
    phylomatching_fits <- infere_anova_phyloanova(
        y = phylomatching_y,
        groups = phylomatching_groups, tree = tree
    )
    random_fits <- infere_anova_phyloanova(
        y = random_y,
        groups = random_groups, tree = tree
    )

    #  pvalues
    phylomatching_all_methods_df <- do.call("rbind", lapply(c("vanilla", "satterthwaite", "lrt"), function(method) {
        phylomatching_pvalues_df <- pvalues_from_fits(
            fit_anova = phylomatching_fits$anova,
            fit_phylolm = phylomatching_fits$phyloanova, tested_method = method,
            tree = tree, REML = FALSE
        )
        phylomatching_pvalues_df$tested_method <- method
        phylomatching_pvalues_df
    }))

    random_all_methods_df <- do.call("rbind", lapply(c("vanilla", "satterthwaite", "lrt"), function(method) {
        random_pvalues_df <- pvalues_from_fits(
            fit_anova = random_fits$anova,
            fit_phylolm = random_fits$phyloanova, tested_method = method,
            tree = tree, REML = FALSE
        )
        random_pvalues_df$tested_method <- method
        random_pvalues_df
    }))

    #  Differientiating the two dataframes before merging
    phylomatching_all_methods_df$group_type <- "phylomatching"
    random_all_methods_df$group_type <- "random"
    data_all_methods_df <- rbind(phylomatching_all_methods_df, random_all_methods_df)

    # Adding the correct_hypothesis column
    data_all_methods_df$correct_hypothesis <- correct_hypothesis


    # Adding the has selected correctly columns
    data_all_methods_df$phylolm_has_selected_correctly <-
        sapply(data_all_methods_df$pvalue_phylolm, function(pvalue) {
            test_selected_correctly(
                correct_hypothesis = correct_hypothesis,
                pvalue, risk_threshold = risk_threshold
            )
        })

    data_all_methods_df$anova_has_selected_correctly <-
        sapply(data_all_methods_df$pvalue_anova, function(pvalue) {
            test_selected_correctly(
                correct_hypothesis = correct_hypothesis,
                pvalue, risk_threshold = risk_threshold
            )
        })

    data_all_methods_df$sim_id <- sim_id
    return(data_all_methods_df)
}

N <- 500
base_values <- c(0, 1)
sigma2_phylo <- 1
sigma2_measure <- 0.1
risk_threshold <- 0.05

dataset_df <- do.call("rbind", lapply(1:N, function(id) {
    rbind(simulate_all_methods(
        sim_id = id,
        correct_hypothesis = "H0",
        base_values = base_values,
        sigma2_phylo = sigma2_phylo,
        sigma2_measure = sigma2_measure,
        risk_threshold = risk_threshold
    ), simulate_all_methods(
        sim_id = id,
        correct_hypothesis = "H1",
        base_values = base_values,
        sigma2_phylo = sigma2_phylo,
        sigma2_measure = sigma2_measure,
        risk_threshold = risk_threshold
    ))
}))

dataset_df %>%
    group_by(tested_method, group_type) %>%
    summarise(
        phylolm_power =
            mean(phylolm_has_selected_correctly[correct_hypothesis == "H1"]),
        phylolm_typeI_error =
            1 - mean(phylolm_has_selected_correctly[correct_hypothesis == "H0"]),
        anova_power =
            mean(anova_has_selected_correctly[correct_hypothesis == "H1"]),
        anova_typeI_error =
            1 - mean(anova_has_selected_correctly[correct_hypothesis == "H0"])
    )

compare_methods <- function(
    N, base_values, risk_threshold, sigma2_phylo,
    sigma2_measure, stoch_process, methods_to_test = c("vanilla", "satterthwaite"), correct_hypothesis = "H1") {
    if (any(!(methods_to_test %in% c("vanilla", "satterthwaite", "lrt")))) {
        stop("Unknown method to test.")
    }
    #  Generating data for each method
    # TODO Utiliser les mêmes données pour les méthodes
    ##  To compute power
    full_power_data <-
        do.call("rbind", lapply(methods_to_test, function(method) {
            data <- simulate_data(
                N = N,
                base_values = base_values,
                risk_threshold = risk_threshold,
                sigma2_phylo = sigma2_phylo,
                sigma2_measure = sigma2_measure,
                test_method = method,
                stoch_process = stoch_process,
                correct_hypothesis = "H1"
            )$data
            #  Adding a column to identify the approximation method
            data$tested_method <- method
            data$metric_type <- "power"
            data
        }))
    ##  To compute type I error
    full_typeI_data <-
        do.call("rbind", lapply(methods_to_test, function(method) {
            data <- simulate_data(
                N = N,
                base_values = base_values,
                risk_threshold = risk_threshold,
                sigma2_phylo = sigma2_phylo,
                sigma2_measure = sigma2_measure,
                test_method = method,
                stoch_process = stoch_process,
                correct_hypothesis = "H0"
            )$data
            #  Adding a column to identify the approximation method
            data$tested_method <- method
            data$metric_type <- "typeI"
            data
        }))

    data <- rbind(full_power_data, full_typeI_data)
    return(
        list(
            data = data,
            sim_parameters = list(
                N = N,
                base_values = base_values,
                risk_threshold = risk_threshold,
                sigma2_phylo = sigma2_phylo,
                sigma2_measure = sigma2_measure,
                stoch_process = stoch_process
            )
        )
    )
}

plot_simulation_data <- function(data, parameters_string, threshold = 0.95) {
    plot_data <- data %>%
        group_by(anova_model, group_type) %>%
        summarize(power = mean(correctly_selected))

    p <- ggplot(plot_data, aes(x = anova_model, y = power, fill = group_type)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_y_continuous(limits = c(0, 1)) +
        labs(
            title = paste0("Power vs Tested Method (", stoch_process, ") | N = ", N, ";", parameters_string),
            x = "Tested Method",
            y = "Power"
        ) +
        geom_hline(yintercept = threshold) +
        theme_minimal()
    p

    return(p)
}

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
