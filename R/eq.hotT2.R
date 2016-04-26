#' Hotelling T2
#' 
#' This function does the one sample Hotelling T2 analog for equivalence tests, determining if a sample's 
#' set of means on repeated measures are practically equivalent, given some pre-specified interval in the metric of Mahalanobis distance.
#' The null hypothesis of nonequality among the means of the responses is rejected if the T2 statistic is less than the critical value.
#' @aliases eq.hotT2
#' @param data a data.frame object
#' @param repeated a character vector containing the names of the repeated measures variables
#' @param ei equivalence interval, in the metric of Mahalanobis distance
#' @param alpha alpha/significance level
#' @references Wellek, S. (2010). \emph{Testing statistical hypotheses of equivalence and noninferiority}. CRC Press.
#' @export eq.hotT2
#' @examples
#' dat <- simRepDat()
#' eq.hotT2(data=dat, repeated=c('repA', 'repB', 'repC'), ei=0.25)


eq.hotT2 <- function(data, repeated, 
    ei, alpha = 0.05) {
    if (any(!complete.cases(data))) {
        print("Missing data present. Only complete cases are used.")
        data <- data[complete.cases(data), 
            ]
    }
    if (class(data) != "data.frame") 
        stop("Data input is not a dataframe.")
    dat <- data[, repeated]
    n <- nrow(dat)
    k <- length(repeated)
    sigma <- cov(dat)
    means <- as.matrix(apply(dat, 2, 
        mean))
    # start hotelling
    check_equiv <- NA
    
    contrasts <- getContrast(k, type = "adjacent")
    vector_of_mean_differences <- contrasts %*% 
        means
    vector_of_mean_differences_pr <- t(vector_of_mean_differences)
    sigma_Mahal <- contrasts %*% sigma %*% 
        t(contrasts)  #covar. matrix with respect to the Mahalanobis distance (vector of mean diffs) 
    sigma_Mahal_inverse <- solve(sigma_Mahal)
    T2 <- n * vector_of_mean_differences_pr %*% 
        sigma_Mahal_inverse %*% vector_of_mean_differences
    fcrit <- (((n - 1) * (k - 1))/(n - 
        k + 1)) * qf(alpha, k - 1, n - 
        k + 1, n * ei^2)
    ifelse(T2 < fcrit, check_equiv <- "T2 is less than the critical value, so the null hypothesis is rejected in favour of equivalence", 
        check_equiv <- "T2 is greater than the critical value, so the null hypothesis is NOT rejected")
    T2_res <- check_equiv
    res <- list(repeatedMeasures = k, 
        means = t(means), ei = ei, T2 = T2, 
        fcrit = fcrit, Decision = check_equiv, 
        sampleSize = n)
    class(res) <- "eq.hotT2"
    return(res)
}
#' @S3method print eq.hotT2
#' @rdname eq.hotT2
#' @method print eq.hotT2
#' @param x object of class \code{eq.hotT2}
print.eq.hotT2 <- function(x, ...) {
    cat("-------Hotelling T2 test for overall equivalence------\n\n")
    cat("There are", x[[1]], "repeated measures for ", 
        n, " participants.", "\n\n")
    cat("The", x[[1]], "means were ")
    cat(x[[2]])
    cat("\n\n")
    cat("The equivalence region was ", 
        x[[3]], "in Mahalanobis distance metric.", 
        "\n\n")
    cat("The T2 statistic was ", x[[4]], 
        "\n\n")
    cat("The F critical value was ", 
        x[[5]], "\n\n")
    cat(x[[6]], "\n")
} 
