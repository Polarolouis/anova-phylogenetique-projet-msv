#!/usr/bin/env Rscript
library(phylolm)
library(phylotools)
library(phytools)
library(ape)
library(ggplot2)
set.seed(1234)

N <- 100 # Number of different simulations
n <- 100

# Arbre
tree <- rphylo(n, 0.1, 0)

## Groupes
K <- 3
get_group <- function(tip) {
    if (tip %in% getDescendants(tree, 105)) {
        return(2)
    }
    if (tip %in% getDescendants(tree, 110)) {
        return(3)
    }
    return(1)
}

source("./simulations/functions-anova.R")

args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
    stop("At least one argument must be supplied", call. = FALSE)
} else if (length(args) >= 3) {
    stoch_process <- args[1]
    sigma2_measure_err <- args[2]
    sigma2_intra_species <- args[3]
    if (length(args)==4) {
        sub_branches <- args[4]
    }
}

# Computing groups
phylo_groups <- as.factor(sapply(1:n, get_group))
non_phylo_groups <- as.factor(sample(c(1, 2, 3), n, replace = TRUE))

calcul_puissance <- function(data, test_method) {
    mean(data[which(data$tested_method == test_method), ]$has_selected_correctly)
}