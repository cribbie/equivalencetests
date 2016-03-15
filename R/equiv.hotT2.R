#' Hotelling T2
#' This function does the one sample Hotelling T2 analog for equivalence tests, determining if a sample's 
#' means on repeated measures are practically equivalent, given some pre-specified interval. 
#' @aliases equiv.hotT2
#' @param data a data.frame object
#' @param repeated a character vector of column names that are the dependent samples
#' @param ei equivalence interval, in the metric of Mahalanobis distance
#' @param alpha alpha level
#' @keywords
#' @export equiv.hotT2
#' @examples
#' dat <- data.frame(pre=rnorm(6), post=rnorm(6,2), fu=rnorm(6,0.5), sex=c(rep('m', 2), rep('f', 4)) )
#' equiv.hotT2()

equiv.hotT2 <- function(data, repeated, ei, alpha = 0.05) {
    if (class(data) != "data.frame") 
        stop("Data input is not a dataframe.")
    dat <- data[, repeated]
    n <- nrow(dat)
    k <- length(repeated)
    sigma <- cov(dat)
    means <- as.matrix(apply(dat, 2, mean))
    # start hotelling
    check_equiv <- NA
    
    contrasts <- getContrast(k, type = "adjacent")
    vector_of_mean_differences <- contrasts %*% 
        means
    vector_of_mean_differences_pr <- t(vector_of_mean_differences)
    sigma_Mahal <- contrasts %*% sigma %*% t(contrasts)  #covar. matrix with respect to the Mahalanobis distance (vector of mean diffs) 
    sigma_Mahal_inverse <- solve(sigma_Mahal)
    T2 <- n * vector_of_mean_differences_pr %*% 
        sigma_Mahal_inverse %*% vector_of_mean_differences
    fcrit <- (((n - 1) * (k - 1))/(n - k + 1)) * 
        qf(alpha, k - 1, n - k + 1, n * ei^2)
    ifelse(T2 < fcrit, check_equiv <- "The null hypothesis is rejected in favour of equivalence", 
        check_equiv <- "The null hypothesis is NOT rejected")
    T2_res <- check_equiv
    res <- list(repeatedMeasures = paste(k, "repeated measures"), 
        means = t(means), ei = paste(ei, "in Mahalanobis distance metric"), 
        T2 = T2, fcrit = fcrit, Decision = check_equiv)
    print(res)
} 
