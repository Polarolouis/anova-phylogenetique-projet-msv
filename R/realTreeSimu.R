# TODO with real tree
# Source : https://genome.cshlp.org/content/29/1/53
# Normaliser l'arbre
require("ape")
require("phylolm")
require("phytools")

tree <- read.tree("./R/chen2019.tree")
tree <- read.tree("./data/data_TER/data_raw/mammals_tree_chen2019.tree")

plotTree(tree)

# Mus et Rat vs le reste

# Groupes équilibrës

