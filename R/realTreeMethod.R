# TODO appliquer Satterthwaite, et calcule les pvalues pour les 5000 genes
# avec Correction Bonferroni et Benjamini-Hochberg (voir Livre Christophe Giraud)

### Import et fonctions utiles
# Repartir du fichier d'analyse Rmd
# Utiliser data.trans, ligne 883 voir RMD
require(phylotools)
require(phytools)
require(phylolm)
require(limma)
require(edgeR)
require(here)
require(ggplot2)
require(dplyr)
require(tidyr)

source("R/utils.R")

### Data import
cdata <- readRDS(here("data", "data_TER", "data", "chen2019_rodents_cpd.rds"))
is.valid <- compcodeR:::check_phyloCompData(cdata)
if (!(is.valid == TRUE)) stop("Not a valid phyloCompData object.")

# Design
design_formula <- as.formula(~condition)
design_data <- compcodeR:::sample.annotations(cdata)[, "condition", drop = FALSE]
design_data$condition <- factor(design_data$condition)
design <- model.matrix(design_formula, design_data)

# Normalisation
nf <- edgeR::calcNormFactors(compcodeR:::count.matrix(cdata) / compcodeR:::length.matrix(cdata), method = "TMM")
lib.size <- colSums(compcodeR:::count.matrix(cdata) / compcodeR:::length.matrix(cdata)) * nf
data.norm <- sweep((compcodeR:::count.matrix(cdata) + 0.5) / compcodeR:::length.matrix(cdata), 2, lib.size + 1, "/")
data.norm <- data.norm * 1e6

# Transformation
data.trans <- log2(data.norm)
rownames(data.trans) <- rownames(compcodeR:::count.matrix(cdata))

### Pvalues computation

#  computing pvalues vec for all genes
pvalue_vec_vanilla <- sapply(seq(1, nrow(data.trans)), function(row_id) {
    trait <- data.trans[row_id, ]
    fit_phylo <- phylolm(trait ~ design_data$condition, phy = cdata@tree, measurement_error = TRUE)
    compute_vanilla_pvalue(fit_phylo)
})

pvalue_vec_vanilla <- setNames(pvalue_vec_vanilla, rownames(data.trans))

pvalue_vec_vanilla_adj <- p.adjust(pvalue_vec_vanilla, method = "BH")

pvalue_vec_satterthwaite <- sapply(seq(1, nrow(data.trans)), function(row_id) {
    trait <- data.trans[row_id, ]
    fit_phylo <- phylolm(trait ~ design_data$condition, phy = cdata@tree, measurement_error = TRUE)
    compute_satterthwaite_pvalue(fit_phylo, tree = cdata@tree)
})

pvalue_vec_satterthwaite <- setNames(pvalue_vec_satterthwaite, rownames(data.trans))

pvalue_vec_satterthwaite_adj <- p.adjust(pvalue_vec_satterthwaite, method = "BH")


pvalue_vec_lrt <- sapply(seq(1, nrow(data.trans)), function(row_id) {
    trait <- data.trans[row_id, ]
    fit_phylo <- phylolm(trait ~ design_data$condition, phy = cdata@tree, measurement_error = TRUE)
    compute_lrt_pvalue(fit_phylo, tree = cdata@tree)
})

pvalue_vec_lrt <- setNames(pvalue_vec_lrt, rownames(data.trans))
pvalue_vec_lrt_adj <- p.adjust(pvalue_vec_lrt, method = "BH")

# REML
pvalue_vec_satterthwaite.REML <- sapply(seq(1, nrow(data.trans)), function(row_id) {
    trait <- data.trans[row_id, ]
    fit_phylo <- phylolm(trait ~ design_data$condition, phy = cdata@tree, measurement_error = TRUE)
    compute_satterthwaite_pvalue(fit_phylo, tree = cdata@tree, REML = TRUE)
})

pvalue_vec_satterthwaite.REML <- setNames(pvalue_vec_satterthwaite.REML, rownames(data.trans))

pvalue_vec_satterthwaite_adj.REML <- p.adjust(pvalue_vec_satterthwaite.REML, method = "BH")

# TODO Histogramme des pvalues

# TODO utiliser UpSetR pour diagramme de Venn

# TODO comparer avec le package evemodel, twothetatest
# Comparer avec OU lrt
# Arbre sans les replicats et les genes data
remotes::install_gitlab("sandve-lab/evemodel")


# computing pvalues vec for all genes
pvalue_vec_vanilla <- sapply(seq(1,nrow(data.trans)), function(row_id) {
    trait <- data.trans[row_id,]
    fit_phylo <- phylolm(trait ~ design_data$condition, phy = cdata@tree, 
    measurement_error = TRUE, model = "OUfixedRoot")
    compute_vanilla_pvalue(fit_phylo)
})

# TODO Utiliser les infos de la ligne 83 du Rmd

# TODO Afficher avec UpSetR les genes differentiellement exprimées et 
# voir les diagrammes de Venn


# Appliquer notre méthode autant de fois que de gène et corriger les pvalues 
# obtenues par la correction pour obtenir

# Vérifier que la F stat = T stat ^ 2