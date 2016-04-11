
#' Internal helper function to compute the gamma Winsorized variance for
# the data in the vector x.  tr is the amount of Winsorization which
# defaults to .2.

winvar <- function(x, tr = 0.2, na.rm = FALSE) {
    if (na.rm) 
        x <- x[!is.na(x)]
    y <- sort(x)
    n <- length(x)
    ibot <- floor(tr * n) + 1
    itop <- length(x) - ibot + 1
    xbot <- y[ibot]
    xtop <- y[itop]
    y <- ifelse(y <= xbot, xbot, y)
    y <- ifelse(y >= xtop, xtop, y)
    winvar <- var(y)
    winvar
}


# create contrast matrices (for hotelling t2: adjacenet mean
# differences; otherwise, pairwise)
getContrast <- function(k, type) {
    if (type == "adjacent") {
        M <- matrix(0, nrow = k - 1, ncol = k)
        for (i in 1:k - 1) {
            M[i, i] <- -1
            M[i, i + 1] <- 1
        }
    }
    if (type == "allPW") {
        M <- matrix(0, nrow = k, ncol = k * (k - 1)/2)
        comb <- combn(k, 2)
        M[cbind(comb[1, ], 1:(k * (k - 1)/2))] <- -1
        M[cbind(comb[2, ], 1:(k * (k - 1)/2))] <- 1
        M <- t(M)
    }
    return(M)
}


#' 'Internal helper Function' for pairwise IUT based tests.
#' 
#' This computes all pairwise mean differences.
#'
#' pairwise_meanDiffs()
pairwise_meanDiffs <- function(sample_means, allcontrasts) {
    mean_diffs <- t(allcontrasts %*% sample_means)  #this is the vector of mean differences, ie. mean difference @ k=1, mean of k=2, etc. Matches mean D, pg 245, table 8.5 
    mean_diff_names <- data.frame(mean_diffs)
    v <- 1:length(sample_means)
    allPairs <- combn(length(v), 2)  # choose a pair from 1:length(v)
    names(mean_diff_names) <- aaply(combn(length(v), 2), 2, function(x) paste0(x[2], 
        "-", x[1]))  # iterate over all pairs
    return(mean_diff_names)
}

#' 'Internal helper Function' for pairwise IUT based tests
#' pairwise_sd()

pairwise_sd <- function(allcontrasts, sigma) {
    
    var_covar <- allcontrasts %*% sigma %*% t(allcontrasts)  #var-covar matrix of intraindiv differences    
    sqrt_varcovar <- sqrt(diag(var_covar))  #SDs of all possible pairwise mean diffs    
    return(sqrt_varcovar)
}

#' 'Internal helper Function' to simulate sample data
#' simRanIntSlope()
simRanIntSlope <- function(sample_size, number_of_sample, equiv_interval, 
    power) {
    id <- rep(1:sample_size, each = number_of_sample)
    gamma00 <- 5  #avg initial status 
    
    # d <- weaken_power_by * equiv_interval #as weaken_power_by increases, d
    # increases.
    
    if (power == FALSE) {
        gamma10 <- equiv_interval  #avg slope      
    } else {
        gamma10 <- equiv_interval  #- d        
    }
    
    timeij <- rep(c(0:(number_of_sample - 1)), times = sample_size)
    
    zeta0i <- rnorm(sample_size, mean = 0, sd = 1)  #int residuals. Int_sd is variance in intercepts. 
    zeta1i <- rnorm(sample_size, mean = 0, sd = 1)  #slope residuals. slope_sd is variance in slopes. 
    eij <- rnorm(sample_size * number_of_sample, mean = 0, sd = 1)  #individual residuals 
    
    y <- gamma00 + gamma10 * timeij + zeta0i + zeta1i * timeij + eij
    # tapply(y,timeij,mean)
    newdat <- data.frame(id, timeij, y)
    return(newdat)
}

#' 'Internal helper Function' to simulate sample data
#' simDat()
simDat <- function() {
    sim_dependent_samples <- function(sample_size, number_of_sample, sigma_matrix, 
        eps, power) {
        if (power == TRUE) {
            population_means <- rep(0, number_of_sample)
        }
        if (power == FALSE) {
            boundary_diff <- eps
            # population_means <- c(0, boundary_diff, rep(0, number_of_sample-2))
            population_means <- c(0, boundary_diff, rep(boundary_diff/2, 
                number_of_sample - 2))
        }
        sample_data <- rmvnorm(sample_size, mean = population_means, sigma = sigma_matrix)
        
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
