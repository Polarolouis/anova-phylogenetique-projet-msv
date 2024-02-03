# Phylocomparison tools
library(phylolm)
library(phylotools)
library(phytools)
library(phylolimma)
library(ape)
library(tidyverse)

# Plotting
library(ggplot2)
library(patchwork)

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

N_simulation_typeI_power <- function(N, base_values, sigma2_phylo, sigma2_measure, risk_threshold = 0.05) {
    df <- do.call("rbind", lapply(1:N, function(id) {
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
}

compute_power_typeI <- function(df) {
    df_plot <- df %>%
        group_by(tested_method, group_type) %>%
        summarise(
            phylolm_power =
                mean(phylolm_has_selected_correctly[correct_hypothesis == "H1"]),
            phylolm_typeIerror =
                1 - mean(phylolm_has_selected_correctly[correct_hypothesis == "H0"]),
            anova_power =
                mean(anova_has_selected_correctly[correct_hypothesis == "H1"]),
            anova_typeIerror =
                1 - mean(anova_has_selected_correctly[correct_hypothesis == "H0"])
        )
    return(df_plot)
}

plot_method_comparison <- function(df_plot, title = "") {
    #  Plot and compare
    anova_plot_typeI <- ggplot(df_plot) +
        aes(x = group_type, y = anova_typeIerror / 3, fill = group_type) +
        ylab("Erreur Type I") +
        xlab("Type de groupe") +
        labs(fill = "Type de groupe") +
        scale_y_continuous(limits = c(0, 1)) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = round(anova_typeIerror, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 0.05)

    anova_plot_power <- ggplot(df_plot) +
        # The /3 is there to account that anova data is repeated 3 times
        aes(x = group_type, y = anova_power / 3, fill = group_type) +
        ylab("Puissance") +
        xlab("Type de groupe") +
        labs(fill = "Type de groupe") +
        scale_y_continuous(limits = c(0, 1)) +
        ggtitle("ANOVA") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = round(anova_power, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9))

    phylolm_plot_typeI <- ggplot(df_plot) +
        aes(x = group_type, y = phylolm_typeIerror, fill = group_type) +
        ylab("Erreur Type I") +
        xlab("Type de groupe") +
        labs(fill = "Type de groupe") +
        scale_y_continuous(limits = c(0, 1)) +
        geom_bar(stat = "identity") +
        geom_text(aes(label = round(phylolm_typeIerror, digits = 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
        geom_hline(yintercept = 0.05) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        facet_wrap(~tested_method)

    phylolm_plot_power <- ggplot(df_plot) +
        aes(x = group_type, y = phylolm_power, fill = group_type) +
        ylab("Puissance") +
        xlab("Type de groupe") +
        labs(fill = "Type de groupe") +
        scale_y_continuous(limits = c(0, 1)) +
        ggtitle("phylolm") +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        facet_wrap(~tested_method) +
        geom_text(aes(label = round(phylolm_power, digits = 2)), vjust = 0.6, position = position_dodge(width = 0.9))

    ((anova_plot_power + phylolm_plot_power + plot_layout(axis_titles = "collect")) / (anova_plot_typeI + phylolm_plot_typeI + plot_layout(axis_titles = "collect"))) + plot_layout(guides = "collect", axes = "collect", axis_titles = "collect") +
        plot_annotation(title = title)
}

## Standardized parameters
total_variance <- 1.0 # sigma2_phylo + sigma2_error, fixed [as tree_height = 1]
heri <- c(0.0, 0.25, 0.5, 1.0) # heritability her = sigma2_phylo / total_variance. 0 means only noise. 1 means only phylo.
snr <- 1 # signal to noise ratio snr = size_effect / total_variance

## Try several parameter values
ggsave <- function(..., bg = "white") ggplot2::ggsave(..., bg = bg)
for (her in heri) {
    sim <- N_simulation_typeI_power(N,
        base_values = c(0, snr * total_variance), 
        sigma2_phylo = her * total_variance,
        sigma2_measure = (1 - her) * total_variance,
    )
    df_sim_plot <- compute_power_typeI(df = sim)

    res_sim_plot <- plot_method_comparison(df_sim_plot, title = paste("BM héritabilité ", her))
    res_sim_plot
    ggsave(paste0("img/simulation_BM_her_", her, ".png"), plot = res_sim_plot)
}
