
simDat <- function() {
  sim_dependent_samples <- function(sample_size, number_of_sample, 
                                    sigma_matrix, eps, power) {
    if (power == TRUE) {
      population_means <- rep(0, number_of_sample)
    }
    if (power == FALSE) {
      boundary_diff <- eps
      # population_means <- c(0, boundary_diff, rep(0,
      # number_of_sample-2))
      population_means <- c(0, boundary_diff, rep(boundary_diff/2, 
                                                  number_of_sample - 2))
    }
    sample_data <- rmvnorm(sample_size, mean = population_means, 
                           sigma = sigma_matrix)
    
    return(sample_data)
  }
  
  make_sigma_matrix <- function(corMatrix, standard_deviation_of_each_sample, 
                                number_of_sample) {
    sigma_matrix <- cor2cov(corMatrix, sd = rep(standard_deviation_of_each_sample, 
                                                number_of_sample))
    return(sigma_matrix)
  }
  k = 3
  times <- 1:k
  rho <- 0.5
  sigma <- 1
  H <- abs(outer(times, times, "-"))
  V <- sigma * rho^H
  p <- nrow(V)
  V[cbind(1:p, 1:p)] <- V[cbind(1:p, 1:p)] * sigma  #cov matrix 
  sigma_matrix <- V  #corMatrix <- cov2cor(V)   
  data <- sim_dependent_samples(sample_size = 100, number_of_sample = k, 
                                sigma_matrix = sigma_matrix, eps = 1, power = TRUE)
  data <- data.frame(data)
  names(data) <- c("repA", "repB", "repC")
  return(data)
} 
