library(ape)
library(phylolm)
library(phytools)

set.seed(4568)

## Tree
n <- 30
tree <- rphylo(n, birth = 0.5, death = 0)

plot(tree, show.tip.label = FALSE, no.margin = TRUE)
nodelabels()

# groups
K <- 3
get_group <- function(tip) {
  if (tip %in% getDescendants(tree, 34)) {
    return(2)
  }
  if (tip %in% getDescendants(tree, 38)) {
    return(3)
  }
  return(1)
}

group <- as.factor(sapply(1:n, get_group))

plot(tree, show.tip.label = FALSE)
tiplabels(bg = group, pch = 21)

# Trait under H0
y <- 2.0 + rTrait(n = 1, phy = tree, model = "BM",
                  parameters = list(acestral.state = 0, sigma2 = 1))


# phylolm fit
fit <- phylolm(y ~ group, phy = tree)
summary(fit)

# Fisher (naive version : uses the inverse of V, while phylolm never computes it.
V <- vcv(tree, model = "BM")
Vinv <- solve(V)
F_stat_naive <- t(fit$fitted.values - mean(y)) %*% Vinv %*% (fit$fitted.values - mean(y))/(K-1)
F_stat_naive <- F_stat_naive / (t(y - fit$fitted.values) %*% Vinv %*% (y - fit$fitted.values)/(n-K))
F_stat_naive

# Fisher (using r squared)
F_stat <- fit$r.squared / (1 - fit$r.squared) * (n - K) / (K - 1)
F_stat

# p-value
p_value <- 1 - pf(F_stat, K - 1, n - K)
p_value

## Check with star tree: phylolm and lm should give the same result
tree <- stree(n)
tree$edge.length <- rep(1.0, nrow(tree$edge))
plot(tree)

# phylo
fit <- phylolm(y ~ group, phy = tree)
F_stat <- fit$r.squared / (1 - fit$r.squared) * (n - K) / (K - 1)
1 - pf(F_stat, K - 1, n - K)

# non phylo
fit <- lm(y ~ group)
aa <- anova(fit)
aa$`F value`
aa$`Pr(>F)`
