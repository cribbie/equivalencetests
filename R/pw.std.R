
#' Pairwise standardized
#' 
#' An equivalence test for repeated measures. Tests the null hypothesis that at least one pairwise mean difference 
#' is not practically equivalent, as defined by an equivalence interval. This null hypothesis is constructed as an intersection-union test, so the null hypothesis fails to be rejected if even one pairwise mean difference is not found to be statistically equivalent.
#' @param data a data.frame object
#' @param repeated a character vector containing the names of the repeated measures variables
#' @param ei equivalence interval, in standardized metric
#' @param alpha  alpha/significance level
#' @export pw.std
#' @examples
#' dat <- simRepDat()
#' pw.std(data=dat, repeated=c('repA', 'repB', 'repC'), ei=1)


pw.std <- function(data, repeated, ei, alpha = 0.05) {
    if (class(data) != "data.frame") 
        stop("Data input is not a dataframe.")
    dat <- data[, repeated]
    n <- nrow(dat)
    k <- length(repeated)
    sigma <- cov(dat)
    means <- as.matrix(apply(dat, 2, mean))
    
    allcontrasts <- getContrast(k, type = "allPW")
    mean_diff_names <- pairwise_meanDiffs(means, 
        allcontrasts)
    sqrt_varcovar <- pairwise_sd(allcontrasts, sigma)
    
    omnibus_res <- NA
    leftside <- NA
    fcrit <- sqrt(qf(alpha, 1, n - 1, n * ei^2))
    for (i in 1:length(mean_diff_names)) {
        leftside <- sqrt(n) * abs((mean_diff_names))/(sqrt_varcovar)
    }
    
    leftside <- unlist(leftside)
    find_nonequiv_res <- which(ifelse(leftside < 
        fcrit, check_equiv <- 1, check_equiv <- 0) == 
        0)
    ifelse(length(find_nonequiv_res) > 0, omnibus_std_res <- "At least one pairwise comparison was not statistically equivalent. There is no evidence for overall equivalence among the repeated measures.", 
        omnibus_std_res <- "All pairwise comparisons were statistically equivalent. There is evidence for overall equivalence among the repeated measures.")  #if at least one pair is NOT equiv, omnibus is not signif. 
    
    if (length(find_nonequiv_res) > 0) {
        print("The following pairwise mean contrasts were found to be nonequivalent:")
        cat("\n")
        print(names(find_nonequiv_res))
        cat("\n\n")
    }
    res <- list(repeatedMeasures = k, means = t(means), 
        ei = ei, Decision = omnibus_std_res)
    class(res) <- "pw.std"
    return(res)
}

#' @rdname pw.std
#' @param x object of class \code{pw.std}
#' @export
print.pw.std <- function(x, ...) {
    cat("------Pairwise unstandardized test for overall equivalence------\n\n")
    cat("There are", x[[1]], "repeated measures.", 
        "\n\n")
    cat("The", x[[1]], "means were ")
    cat(x[[2]])
    cat("\n\n")
    cat("The equivalence interval was ", x[[3]], 
        "in standardized metric.", "\n\n")
    print(x[[4]])
} 
