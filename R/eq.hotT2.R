#' Hotelling T2
#' 
#' This function does the one sample Hotelling T2 analog for equivalence tests, determining if a sample's 
#' set of means on repeated measures are practically equivalent, given some pre-specified interval in the metric of Mahalanobis distance.
#' @aliases eq.hotT2
#' @param data a data.frame object
#' @param repeated a character vector containing the names of the repeated measures variables
#' @param ei equivalence interval, in the metric of Mahalanobis distance
#' @param alpha alpha/significance level
#' @references Wellek, S. (2010). \emph{Testing statistical hypotheses of equivalence and noninferiority}. CRC Press.
#' @export eq.hotT2
#' @examples
#' k.obs <- 3
#' elements <- c(1, 0.7, 0.8, 1, 0.5,1)  
#' 
#' X <- diag(k.obs)  
#' X[lower.tri(X, diag = TRUE)] <- elements
#' X <- X + t(X) - diag(diag(X)) 
#' colnames(X) <- c('repA', 'repB', 'repC')
#' rownames(X) <- c('repA', 'repB', 'repC')
#' 
#' if(!require(mvtnorm)) install.packages(mvtnorm); library(mvtnorm) 
#' dat <- data.frame(rmvnorm(n=40, sigma=X))
#' names(dat) <- c('repA', 'repB', 'repC')
#' 
#' eq.hotT2(data=dat, repeated=c('repA', 'repB', 'repC'), ei=0.25)


eq.hotT2 <- function(data, repeated, ei, alpha = 0.05) {
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
    vector_of_mean_differences <- contrasts %*% means
    vector_of_mean_differences_pr <- t(vector_of_mean_differences)
    sigma_Mahal <- contrasts %*% sigma %*% t(contrasts)  #covar. matrix with respect to the Mahalanobis distance (vector of mean diffs) 
    sigma_Mahal_inverse <- solve(sigma_Mahal)
    T2 <- n * vector_of_mean_differences_pr %*% sigma_Mahal_inverse %*% 
        vector_of_mean_differences
    fcrit <- (((n - 1) * (k - 1))/(n - k + 1)) * qf(alpha, k - 
        1, n - k + 1, n * ei^2)
    ifelse(T2 < fcrit, check_equiv <- "The null hypothesis is rejected in favour of equivalence", 
        check_equiv <- "The null hypothesis is NOT rejected")
    T2_res <- check_equiv
    res <- list(repeatedMeasures = paste(k, "repeated measures"), 
        means = t(means), ei = paste(ei, "in Mahalanobis distance metric"), 
        T2 = T2, fcrit = fcrit, Decision = check_equiv)
    print(res)
} 
