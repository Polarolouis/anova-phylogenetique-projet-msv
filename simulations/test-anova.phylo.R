library(phylolm)
library(phylotools)
library(phytools)
library(ape)
set.seed(1234)

n <- 100
tree <- rphylo(n, 0.1, 0)

sigma2err <- 1

# Continuous phylo trait
trait <- rTrait(1, tree, "BM")

# Adding noise to the trait
trait <- trait + rnorm(n, mean = 0, sqrt(sigma2err))

# Simulation positive

## Groupes
get_group <- function(tip) {
    if (tip %in% getDescendants(tree, 105)) {
        return(2)
    }
    if (tip %in% getDescendants(tree, 110)) {
        return(3)
    }
    return(1)
}

group <- as.factor(sapply(1:n, get_group))

## Réponse
mu1 <- 2
mu2 <- -5
mu3 <- 2
y <- mu1 * (group == 1) + mu2 * (group == 2) + mu3 * (group == 3)
y <- y + trait

par(mar = c(5, 0, 0, 0) + 0.1)
plot(tree, show.tip.label = FALSE, x.lim = 50)
tiplabels(bg = group, pch = 21)
phydataplot(y, tree, scaling = 0.1, offset = 4)

fitANOVA <- lm(y ~ group)

fitphyloANOVA <- phylolm(y ~ group, phy = tree)

# Simulation négative

groups_non_phylo <- as.factor(sample(c(1,2,3), n, replace = TRUE))
y_non_phy <- mu1 * (groups_non_phylo == 1) + mu2 * (groups_non_phylo == 2) + mu3 * (groups_non_phylo == 3)
y_non_phy <- y_non_phy + trait

par(mar = c(5, 0, 0, 0) + 0.1)
plot(tree, show.tip.label = FALSE, x.lim = 50)
tiplabels(bg = groups_non_phylo, pch = 21)
phydataplot(y_non_phy, tree, scaling = 0.1, offset = 4)

fit_nonphy_ANOVA <- lm(y_non_phy ~ groups_non_phylo)
fitphy_nonphy_ANOVA <- phylolm(y_non_phy ~ groups_non_phylo, phy = tree)

# Summary
summary(fit_nonphy_ANOVA)
summary(fitphy_nonphy_ANOVA)