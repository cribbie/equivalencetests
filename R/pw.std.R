
#' Pairwise standardized
#' 
#' An equivalence test for repeated measures. Tests the null hypothesis that at least one pairwise mean difference 
#' is not practically equivalent, as defined by an equivalent interval.  
#' @param data data set in data.frame format
#' @param repeated a character vector of column names designating the repeated measures
#' @param ei equivalence interval 
#' @param alpha  alpha level
#' @export pw.std
#' @examples
#' pw.std()


pw.std <- function(data, repeated, ei, alpha = 0.05) {
    if (class(data) != "data.frame") 
        stop("Data input is not a dataframe.")
    dat <- data[, repeated]
    n <- nrow(dat)
    k <- length(repeated)
    sigma <- cov(dat)
    means <- as.matrix(apply(dat, 2, mean))
    
    allcontrasts <- getContrast(k, type = "allPW")
    mean_diff_names <- pairwise_meanDiffs(means, allcontrasts)
    sqrt_varcovar <- pairwise_sd(allcontrasts, sigma)
    
    omnibus_res <- NA
    leftside <- NA
    fcrit <- sqrt(qf(alpha, 1, n - 1, n * ei^2))
    for (i in 1:length(mean_diff_names)) {
        leftside <- sqrt(n) * abs((mean_diff_names))/(sqrt_varcovar)
    }
    
    leftside <- unlist(leftside)
    find_nonequiv_res <- which(ifelse(leftside < fcrit, check_equiv <- 1, 
        check_equiv <- 0) == 0)
    ifelse(length(find_nonequiv_res) > 0, omnibus_std_res <- "No evidence for equivalence", 
        omnibus_std_res <- "Evidence for equivalence")  #if at least one pair is NOT equiv, omnibus is not signif. 
    
    res <- list(repeatedMeasures = paste(k, "repeated measures"), 
        means = t(means), ei = paste(ei, "in standardized metric"), 
        Decision = omnibus_std_res)
    print(res)
    return(res)
} 
