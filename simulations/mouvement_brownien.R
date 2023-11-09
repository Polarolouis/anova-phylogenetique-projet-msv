library(ggplot2)

# Number of BM indep coordinates
N <- 1000
time_delta <- 0.5

generate_brownian_movement_traj <- function(N, time_delta) {
    # BM indep vector
    bm_indep <- rnorm(N, mean = 0, time_delta)

    # The trajectory
    bm_traj <- cumsum(bm_indep)
    time_serie <- cumsum(rep(time_delta, N))

    return(data.frame(time_serie = time_serie, bm_traj = bm_traj))
}

bm_traj_and_time <- generate_brownian_movement_traj(N, time_delta)

# Plotting
ggplot(bm_traj_and_time) +
    aes(x = time_serie, y = bm_traj) +
    geom_line()
