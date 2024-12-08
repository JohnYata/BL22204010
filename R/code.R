
#' @importFrom Rcpp sourceCpp

#' @importFrom Rcpp evalCpp
#' Simulate a Multi-Armed Bandit Environment
#'
#' This function simulates a multi-armed bandit environment with Gaussian rewards and non-stationary means.
#'
#' @param N Integer. The number of arms in the bandit.
#' @param T Integer. The number of time steps for the simulation.
#' @param sigma Numeric. The standard deviation of the rewards.
#' @param delta Numeric. The standard deviation of the mean updates.
#' @param mu Numeric vector. The initial means for each arm.
#' @return A numeric matrix containing the simulated rewards for each time step and arm.
#' @export
#' @importFrom stats rnorm runif

simulate_bandit <- function(N, T, sigma, delta) {

  mu <- runif(N, min = 0, max = 1)
  rewards <- matrix(0, nrow = T, ncol = N)
  means <- matrix(0, nrow = T, ncol = N)
  means[1, ] <- mu

  for (t in 1:T) {
    for (i in 1:N) {

      rewards[t, i] <- rnorm(1, mean = mu[i], sd = sigma)
      if (t < T) {
        mu[i] <- mu[i] + rnorm(1, mean = 0, sd = delta)
        means[t + 1, i] <- mu[i]
      }
    }
  }
  list(r = rewards, means = means)
}


