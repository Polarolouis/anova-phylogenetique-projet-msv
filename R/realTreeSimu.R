# TODO with real tree
# Source : https://genome.cshlp.org/content/29/1/53
# Normaliser l'arbre
require("ape")
require("phylolm")
require("phytools")

source("./R/utils.R")

plot_group_on_tree <- function(tree, groups) {
    plot(tree)
    tiplabels(bg = groups, pch = 21)
}

tree <- read.tree("./R/chen2019.tree")
# Normalising tree edge length
taille_tree <- diag(vcv(tree))[1]
tree$edge.length <- tree$edge.length / taille_tree

plotTree(tree, ftype="i")

# Mus et Rat vs le reste

group_mus_rat_vs_other <- sapply(44:(44+41), function(tip) {
    if (tip %in% getDescendants(tree = tree, 55)) {
        return(1)
    }
    return(2)
})

plot_group_on_tree(tree, group = group_mus_rat_vs_other)


# Groupes équilibrës
