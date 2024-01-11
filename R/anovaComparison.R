# Phylocomparison tools
library(phylolm)
library(phylotools)
library(phytools)
library(phylolimma)
library(ape)

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

# Compute traits
# Parameters
base_values <- c(1, 2.1) # The base trait to add
sigma2_phylo <- 1
sigma2_measure <- 0
stoch_process <- "BM"
hypo_test <- "satterthwaite" #"vanilla" # "satterthwaite", "likelihood_ratio"


matching_phylo_traits <- compute_trait_values(
    groups = phylo_matching_groups,
    base_values = base_values, tree,
    sigma2_phylo = sigma2_phylo, sigma2_measure = sigma2_measure,
    stoch_process = stoch_process
)

random_traits <- compute_trait_values(
    groups = random_groups,
    base_values = base_values, tree,
    sigma2_phylo = sigma2_phylo, sigma2_measure = sigma2_measure,
    stoch_process = stoch_process
)

# For phylo matching
matching_anova <- lm(matching_phylo_traits ~ phylo_matching_groups)

# TODO Handling the stoch process and model for phylolm
model <- stoch_process

matching_phyloanova <- phylolm(matching_phylo_traits ~ phylo_matching_groups,
    phy = tree,
    model = model,
    measurement_error = TRUE # To let phylolm know if there's measurement error
)

matching_phyloanova$REML = FALSE

if (hypo_test %in% c("vanilla", "satterthwaite")) {
    phyloanova_F_stat <- compute_F_statistic(
        r_squared = matching_phyloanova$r.squared,
        df1 = K - 1,
        df2 = nb_species - K
    )

    df1 <- K - 1
    df2 <- nb_species - K

    if (hypo_test == "satterthwaite") {
        df2 <- phylolimma:::ddf_satterthwaite(matching_phyloanova, tree)
    }

    p_value <- 1 - pf(phyloanova_F_stat, df1 = df1, df2 = df2)
}

if (hypo_test == "likelihood_ratio") {
    # How to obtain the loglikehood under H0 ?
    # TODO Find the correct way to do this
    # I assume that under H0 this is like saying everyone is from the same group
    h0_matching_phyloanova <- phylolm(matching_phylo_traits ~ rep(1, nb_species),
        phy = tree,
        model = model,
        measurement_error = (sigma2_measure != 0) # To let phylolm know if there's measurement error
    )
    # But this gives a LAPACK error, the system is not inversible.

    lambda_ratio_stat <- -2(h0_matching_phyloanova$logLik - matching_phyloanova$logLik)
}

p_value