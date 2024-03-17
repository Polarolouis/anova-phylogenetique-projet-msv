library(ggplot2)
library(dplyr)

# # Number of BM indep coordinates
# N <- 1000
# time_delta <- 0.5

generate_brownian_motion_traj <- function(N, time_delta, variance = 1) {
    # BM indep vector
    bm_indep <- rnorm(N, mean = 0, sqrt(time_delta * variance))

    # The trajectory
    bm_traj <- cumsum(bm_indep)
    time_serie <- cumsum(rep(time_delta, N))

    return(data.frame(time_serie = time_serie, bm_traj = bm_traj))
}

# bm_traj_and_time <- generate_brownian_motion_traj(N, time_delta, 1)

# # Plotting
# ggplot(bm_traj_and_time) +
#     aes(x = time_serie, y = bm_traj) +
#     geom_hline(aes(yintercept = 0, color = "red")) +
#     geom_line()

# # Generate multiple BM
# n <- 100
# multiple_BM <- do.call(rbind, lapply(seq(1, n, 1), function(idx) {
#     data <- generate_brownian_motion_traj(N, time_delta)
#     data$id <- rep(idx, N)
#     data <- data[, c(3, 1:2)]
#     return(data)
# }))

# ## Plotting multiple BMs
# ggplot(multiple_BM) +
#     aes(x = time_serie, y = bm_traj, group = as.factor(id), alpha = 0.5) +
#     geom_hline(aes(yintercept = 0, color = "red")) +
#     geom_line()

# For phylogenic tree

generate_phylo_tree <- function(n_tips, N, time_delta) {
    # Les instants
    time_serie <- cumsum(rep(time_delta, N))

    # Instants de spéciations
    spec_times <- sort(sample(seq(1:N - 10), n_tips, replace = TRUE))

    # Choisir la branche dont on se sépare, au premier choix 1 seule, au second 2
    branche_orig_spec <- sapply(
        seq_len(length(spec_times)),
        function(nb_choix) sample(seq(1, nb_choix - 1), size = 1)
    )

    # On ajoute alors la première espèce
    spec_times <- c(0, spec_times)

    # Les trajectoires
    liste_spec_traj <- lapply(seq_len(length(spec_times)), function(idx) {
        generate_brownian_motion_traj(
            N = N - spec_times[idx],
            time_delta = time_delta
        )[["bm_traj"]]
    })

    phylo_tree_data <- data.frame(
        id_branche = factor(),
        time = numeric(), traj = numeric(), is_spec_time = logical()
    )

    # Pour la construction récursive de l'arbre
    for (idx in seq_len(length(liste_spec_traj)-1)) {
        # Les k pas manquants
        k <- spec_times[idx]

        if (k == 0) {
            current_data <- data.frame(
                id_branche = rep(
                    idx,
                    length(liste_spec_traj[[idx]])
                ),
                time = time_serie, traj = liste_spec_traj[[idx]], is_spec_time = rep(FALSE, N)
            )
        } else {
            base <- phylo_tree_data[phylo_tree_data[["id_branche"]] == (branche_orig_spec[[idx]]),"traj"][1:k]
            if(anyNA(base)){
                print(spec_times[[idx]])
                print(phylo_tree_data$id_branche)
                print(branche_orig_spec[[idx]])
            }
            traj <-
                # On prend la dernière valeur avant div qu'on ajoute à la 
                # trajectoire
                c(base, tail(base, n = 1) + liste_spec_traj[[idx]])
            current_data <- data.frame(
                id_branche = rep(idx, N),
                time = time_serie, traj = traj, is_spec_time = (time_serie == time_serie[k])
            )
        }

        phylo_tree_data <- rbind(phylo_tree_data, current_data)
    }

    return(phylo_tree_data)
}

# debug(generate_phylo_tree)
