library(ggplot2)

# Number of BM indep coordinates
N <- 1000
time_delta <- 0.5

generate_brownian_motion_traj <- function(N, time_delta, variance = 1) {
    # BM indep vector
    bm_indep <- rnorm(N, mean = 0, sqrt(time_delta * variance))

    # The trajectory
    bm_traj <- cumsum(bm_indep)
    time_serie <- cumsum(rep(time_delta, N))

    return(data.frame(time_serie = time_serie, bm_traj = bm_traj))
}

bm_traj_and_time <- generate_brownian_motion_traj(N, time_delta, 1)

# Plotting
ggplot(bm_traj_and_time) +
    aes(x = time_serie, y = bm_traj) +
    geom_hline(aes(yintercept = 0, color = "red")) +
    geom_line()

# Generate multiple BM
n <- 100
multiple_BM <- do.call(rbind, lapply(seq(1,n,1), function(idx) {
    data <- generate_brownian_motion_traj(N, time_delta)
    data$id <- rep(idx, N)
    data <- data[,c(3,1:2)]
    return(data)
}))

## Plotting multiple BMs
ggplot(multiple_BM) +
    aes(x = time_serie, y = bm_traj, group = as.factor(id), alpha = 0.5) +
    geom_hline(aes(yintercept = 0, color = "red")) +
    geom_line()

# For phylogenic tree

generate_phylo_tree <- function(n_tips, max_time) {
    
}