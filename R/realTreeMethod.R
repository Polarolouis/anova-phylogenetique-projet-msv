# TODO appliquer Satterthwaite, et calcule les pvalues pour les 5000 genes
# avec Correction Bonferroni et Benjamini-Hochberg (voir Livre Christophe Giraud)

### Import et fonctions utiles
# Repartir du fichier d'analyse Rmd
# Utiliser data.trans, ligne 883 voir RMD
library(phylotools)
library(phytools)
library(phylolm)
library(limma)
library(edgeR)
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)

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
    compute_satterthwaite_pvalue(fit_phylo, tree = cdata@tree, )
})

# TODO Analyser l'origine de la surestimation du nombre de df2 pour le gène 1. Vient pê de sigma2_error ~ 1e-11

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
## Préparation du dataframe
pvalues_dataframe <- data.frame(
    gene = rep(rownames(data.trans), 8),
    pvalue = c(pvalue_vec_vanilla, pvalue_vec_vanilla_adj, 
        pvalue_vec_satterthwaite, pvalue_vec_satterthwaite_adj, 
        pvalue_vec_lrt, pvalue_vec_lrt_adj, pvalue_vec_satterthwaite.REML, 
        pvalue_vec_satterthwaite_adj.REML),
    test_method = rep(c("Vanilla", "VanillaAdj", "Satterthwaite", 
    "SatterthwaiteAdj", "LRT", "LRTAdj", "SatterthwaiteREML", 
    "SatterthwaiteAdjREML"), each = nrow(data.trans))
)
pvalues_dataframe$test_method <- as.factor(pvalues_dataframe$test_method)
pvalues_dataframe <- pvalues_dataframe %>% mutate(selected = ifelse(pvalue < 0.05, 1, 0))

pvalues_dataframe_wide <- pvalues_dataframe %>% 
    pivot_wider(id_cols = gene, 
    names_from = test_method, 
    values_from = selected) %>%
    data.frame()

## Graphiques
ggplot(pvalues_dataframe) +
    aes(x = genes, y = pvalues, fill = test_method) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~test_method)

# DONE utiliser UpSetR pour diagramme de Venn
library(UpSetR)
upset(pvalues_dataframe_wide, 
    nsets = 8,
    mainbar.y.label = "Nombre de gènes en commun",
    sets.x.label = "Nombre de gènes sélectionnés")

# TODO comparer avec le package evemodel, twothetatest
# Comparer avec OU lrt
# Arbre sans les replicats et les genes data
remotes::install_gitlab("sandve-lab/evemodel")



# TODO Utiliser les infos de la ligne 83 du Rmd

# TODO Afficher avec UpSetR les genes differentiellement exprimées et
# voir les diagrammes de Venn


# Appliquer notre méthode autant de fois que de gène et corriger les pvalues
# obtenues par la correction pour obtenir

# Vérifier que la F stat = T stat ^ 2