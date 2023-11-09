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

bm_traj_and_time <- generate_brownian_motion_traj(N, time_delta)

# Generate multiple BM


# For phylogenic tree

generate_phylo_tree <- function(n_tips, max_time) {

}