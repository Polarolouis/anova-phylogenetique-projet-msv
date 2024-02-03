# Copyright (C) {2023} {PB, MG, AC}
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

# Program to download and format the data from:
# Chen, J., Swofford, R., Johnson, J., Cummings, B. B., Rogel, N., Lindblad-Toh, 
# K., Haerty, W., Palma, F. di, & Regev, A. (2019). A quantitative framework for 
# characterizing the evolutionary history of mammalian gene expression. 
# Genome Research, 29(1), 53–63. https://doi.org/10.1101/gr.237636.118

## Install correct version of compcodeR
# devtools::install_github("csoneson/compcodeR", ref = "phylocomp")

library(here)
library(compcodeR)

################################################################################
## Read and format count and length data
################################################################################
# Read the raw counts
raw_counts_full <- read.table(here("data/data_TER/data/rawCounts.txt"))
raw_counts <- raw_counts_full[,-c(1,2)]
raw_counts <- round(raw_counts) # expected RSEM values -> counts
## format sample names
vectrep <- paste0(".",sapply(1:dim(raw_counts)[2],function(i) length(strsplit(colnames(raw_counts),"[.]")[[i]]))-1)
colnames(raw_counts) <- paste0(colnames(raw_counts),replace(vectrep,vectrep==".1",""))


# Read length information
leng_full <- read.table(here("data/data_TER/data/rawLengths.txt"))
leng <- leng_full[,-c(1,2)]
## format sample names
vectrep <- paste0(".",sapply(1:dim(leng)[2],function(i) length(strsplit(colnames(leng),"[.]")[[i]]))-1)
colnames(leng) <- paste0(colnames(leng),replace(vectrep,vectrep==".1",""))

# remove NA
raw_counts_noNA <- raw_counts[complete.cases(raw_counts),]
length_noNA <- leng[complete.cases(raw_counts),colnames(raw_counts_noNA)]

# remove zero counts
zero_counts <- rowSums(raw_counts_noNA) == 0
raw_counts_noNA <- raw_counts_noNA[!zero_counts, ]
length_noNA <- length_noNA[!zero_counts, ]

# # remove counts when at least one zero
# zero_counts <- apply(raw_counts_noNA, 1, function(x) any(x == 0))
# raw_counts_noNA <- raw_counts_noNA[!zero_counts, ]
# length_noNA <- length_noNA[!zero_counts, ]

################################################################################
## Read and format condition data
################################################################################

mammalsconds <- read.csv2(here("data/data_TER/data_raw/speciesgroups.txt"), sep = "\t")
mammalsconds
## "primates" : 2= primates, 1= non primates
## "rodents" : 2 = rodents , 1 = non rodents

## for primates only... 
condName <- "primates"
# Species condition (1 = non primates, 2 = primates)
condstp <- mammalsconds$primates
names(condstp) <- mammalsconds$species
# Sample Condition
idSpe <- sapply(1:dim(raw_counts)[2],function(i) strsplit(colnames(raw_counts)[i],"[.]")[[1]][1])
primates <- condstp[idSpe]  
names(primates) <- colnames(raw_counts)
# Format
primates <- factor(primates)
colData <- data.frame(primates)
colData$primates <- factor(colData$primates)
colData$species <- sapply(1:length(rownames(colData)),function(i) strsplit(rownames(colData),"[.]")[[i]][1])


## for rodents only... (plus favorable pour philoebayes)
condName <- "rodents"
# Species condition (1 = non rodents, 2 = rodents)
condstp <- mammalsconds$rodents
names(condstp) <- mammalsconds$species
# Sample Condition
rodents <- condstp[idSpe]  
names(rodents) <- colnames(raw_counts)
# Format
colData$rodents <- factor(rodents)

################################################################################
## Read and format Tree ## COPIE COLLE DE STERN CA MARCHE PO 
################################################################################
library(ape)
## Get the tree
tree <- read.tree(file = here("data/data_TER/data_raw/mammals_tree_trimmed_hassler2020.tree"))
tree$node.label <- NULL
tree$tip.label
plot(tree)

## Species names
# Match
tree_data_cor <- match(tree$tip.label, colData$species)
data_tree_cor <- match(colData$species, tree$tip.label)
# Species in the tree NOT in data
tree$tip.label[is.na(tree_data_cor)]
# Species in data NOT in the tree
colData$species[is.na(data_tree_cor)]

## Format Tree
# Get rid of species not in data
tree <- drop.tip(tree, tip = tree$tip.label[is.na(tree_data_cor)])
plot(tree)
# function to add replicates
#createTreeWithReplicates <- function(tree){
  tree_rep <- tree
  for (tip_label in tree$tip.label) {
    print(tip_label)
    replis <- colnames(raw_counts_noNA)[grepl(tip_label, colnames(raw_counts_noNA))]
    print(replis)
    for (rep in replis) {
      print(rep)
      tree_rep <- phytools::bind.tip(tree_rep, tip.label = rep,
                                     where = which(tree_rep$tip.label == tip_label))
      plot(tree_rep)
    }
  }
  # Remove original tips
#  return(ape::drop.tip(tree_rep, tree$tip.label))
#}
tree_rep <- ape::drop.tip(tree_rep, tree$tip.label)
#tree_rep <- createTreeWithReplicates(tree)
# Plot
plot(tree_rep)

## Reorder data
corr <- match(tree_rep$tip.label, rownames(colData))
colData <- colData[corr, ]
raw_counts_noNA <- raw_counts_noNA[, corr]
length_noNA <- length_noNA[, corr]

## Remove if zero count for all replicates of a species
zero_counts <- apply(raw_counts_noNA, 1, function(x) any(tapply(t(x), colData$species, FUN = sum) == 0))
raw_counts_noNA <- raw_counts_noNA[!zero_counts, ]
length_noNA <- length_noNA[!zero_counts, ]


# # ################################################################################
# # ## Read and format Tree ## MODIF ANAELLE, MARCHE PO POUR CompcodeR mais j'ai pas compris pk
# # ################################################################################
# library(ape)
# ## Get the tree
# tree <- read.tree(file = here("data_raw/mammals_tree_trimmed_hassler2020.tree"))
# tree$node.label <- NULL
# tree$tip.label
# plot(tree)
# 
# ## Arbre initial
# col <- colnames(raw_counts)
# col <- sub("[.][[:digit:]]+", "", col)
# 
# tree_rep <- tree # Création de l'arbre avec les réplicats
# 
# for (i in 1:length(tree$tip.label)){
#   print(i)
#   label <- tree$tip.label[i]
#   print(label)
#   nb_replicates <- length(which(col == label))
#   print(nb_replicates)
#   if(nb_replicates>1){               # Pour éviter de rentrer ds la boucle si pas besoin (et éviter une erreur)
#     for (j in (nb_replicates-1):1){  # J'ai inversé l'ordre d'ajout des réplicats pour qu'ils soient ds le bon ordre ds l'arbre
#       new_replicate = paste(label,j,sep = ".")
#       tree_rep <- phytools::bind.tip(tree_rep,new_replicate,edge.length = NULL, which(tree_rep$tip.label==label),position=0)
#     }
#   }
# }
# plot(tree_rep)
# 
# 
# ## Reorder data
# corr <- match(tree_rep$tip.label, rownames(colData))
# colData <- colData[corr, ]
# raw_counts_noNA <- raw_counts_noNA[, corr]
# length_noNA <- length_noNA[, corr]
# 
# 


################################################################################
# Format for compcodeR
################################################################################
## rename key column "condition" in colData
colnames(colData)[grep(condName, colnames(colData))] <- "condition"
## id species 
colData$id.species <- colData$species
## create object
cpd <- phyloCompData(tree = tree_rep,
                     count.matrix = raw_counts_noNA, 
                     sample.annotations = colData, 
                     info.parameters = list(dataset = "chen2019_rodents_cpd", uID = "1"),
                     length.matrix = as.matrix(length_noNA))
## check obj conformity
check_phyloCompData(cpd)

## save data to rds file
saveRDS(cpd, file = here("data/data_TER/data", "chen2019_rodents_cpd.rds"))
