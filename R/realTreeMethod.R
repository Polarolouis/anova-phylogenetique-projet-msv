# TODO appliquer Satterthwaite, et calcule les pvalues pour les 5000 genes
# avec Correction Bonferroni et Benjamini-Hochberg (voir Livre Christophe Giraud)

# Repartir du fichier d'analyse Rmd
# Utiliser data.trans, ligne 883 voir RMD

require(limma)
require(edgeR)
library(here)

cdata <- readRDS(here("data","data_TER","data", "chen2019_rodents_cpd.rds"))
is.valid <- compcodeR:::check_phyloCompData(cdata)
if (!(is.valid == TRUE)) stop('Not a valid phyloCompData object.')

# Design
design_formula <- as.formula(~ condition)
design_data <- compcodeR:::sample.annotations(cdata)[, 'condition', drop = FALSE]
design_data$condition <- factor(design_data$condition)
design <- model.matrix(design_formula, design_data)

# Normalisation
nf <- edgeR::calcNormFactors(compcodeR:::count.matrix(cdata) / compcodeR:::length.matrix(cdata), method = 'TMM')
lib.size <- colSums(compcodeR:::count.matrix(cdata) / compcodeR:::length.matrix(cdata)) * nf
data.norm <- sweep((compcodeR:::count.matrix(cdata) + 0.5) / compcodeR:::length.matrix(cdata), 2, lib.size + 1, '/')
data.norm <- data.norm * 1e6

# Transformation
data.trans <- log2(data.norm)
rownames(data.trans) <- rownames(compcodeR:::count.matrix(cdata))

# computing pvalues vec for all genes
pvalue_vec_vanilla <- sapply(seq(1,nrow(data.trans)), function(row_id) {
    trait <- data.trans[row_id,]
    fit_phylo <- phylolm(trait ~ design_data$condition, phy = cdata@tree)
    compute_vanilla_pvalue(fit_phylo)
})

pvalue_vec_vanilla <- setNames(pvalue_vec_vanilla, rownames(data.trans))

pvalue_vec_vanilla_adj <- p.adjust(pvalue_vec_vanilla, method = "BH")

pvalue_vec_satterthwaite <- sapply(seq(1,nrow(data.trans)), function(row_id) {
    trait <- data.trans[row_id,]
    fit_phylo <- phylolm(trait ~ design_data$condition, phy = cdata@tree)
    compute_satterthwaite_pvalue(fit_phylo, tree = cdata@tree)
})

pvalue_vec_satterthwaite <- setNames(pvalue_vec_satterthwaite, rownames(data.trans))

pvalue_vec_satterthwaite_adj <- p.adjust(pvalue_vec_satterthwaite, method = "BH")


pvalue_vec_lrt <- sapply(seq(1,nrow(data.trans)), function(row_id) {
    trait <- data.trans[row_id,]
    fit_phylo <- phylolm(trait ~ design_data$condition, phy = cdata@tree)
    compute_lrt_pvalue(fit_phylo, tree = cdata@tree)
})

pvalue_vec_lrt <- setNames(pvalue_vec_lrt, rownames(data.trans))
pvalue_vec_lrt_adj <- p.adjust(pvalue_vec_lrt, method = "BH")



# Appliquer notre méthode autant de fois que de gène et corriger les pvalues 
# obtenues par la correction pour obtenir

# Vérifier que la F stat = T stat ^ 2