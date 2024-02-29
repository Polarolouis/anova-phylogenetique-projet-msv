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

# TODO phylolm(stoch_model = "OUfixedRoot")

fit_phylo <- phylolm(trait ~ design_data$condition, phy = cdata@tree, measurement_error = TRUE)
compute_satterthwaite_pvalue(fit_phylo, tree = cdata@tree, return_df = TRUE)
