#' Schuirmann's Test of the Equivalence of Two Independent Groups
#' 
#' Perform Schuirmann's two-sample equivalence test using Yuen's formula for
#' trimmed means on the data in x and y. Missing values (values stored as NA) are 
#' automatically removed.
#' 
#' @aliases yuen_equiv
#' @param x a numeric vector
#' @param y a numeric vector
#' @param equivint equivalence interval
#' @param tr proportion of data to trim
#' @param alpha desired alpha level
#' @param ... additional arguments to be passed
#' 
#' @author Rob Cribbie \email{cribbie@@yorku.ca} and 
#'   Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export yuen_equiv
#' @examples
#' \dontrun{
#' #equivalence correlation test between v1 and v2 with an interval of .2
#' v1 <- rnorm(100)
#' v2 <- v1 + rnorm(100, 2)
#' yuen_equiv(v1, v2, .2)
#' }
equiv_t <- function(x, y, equivint, alpha = 0.05, varequiv = FALSE, na.rm = TRUE, ...) {
    if (na.rm) 
        x <- x[!is.na(x)]
    if (na.rm) 
        y <- y[!is.na(y)]
    if (varequiv == FALSE) {
        t1 <- (mean(x) - mean(y) - equivint)/sqrt((var(x)/length(x)) + (var(y)/length(y)))
        t2 <- (mean(x) - mean(y) + equivint)/sqrt((var(x)/length(x)) + (var(y)/length(y)))
        dft <- (((var(x)/length(x)) + (var(y)/length(y)))^2)/((var(x)^2/(length(x)^2 * (length(x) - 1))) + 
            (var(y)^2/(length(y)^2 * (length(y) - 1))))
        probt1 <- pt(t1, dft, lower.tail = T)
        probt2 <- pt(t2, dft, lower.tail = F)
        ifelse(probt1 < alpha & probt2 < alpha, decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected", 
            decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected")
    }
    if (varequiv == TRUE) {
        t1 <- (mean(x) - mean(y) - equivint)/sqrt(((((length(x) - 1) * sd(x)^2) + ((length(y) - 1) * sd(y)^2))/(length(x) + 
            length(y) - 2)) * (1/length(x) + 1/length(y)))
        t2 <- (mean(x) - mean(y) + equivint)/sqrt(((((length(x) - 1) * sd(x)^2) + ((length(y) - 1) * sd(y)^2))/(length(x) + 
            length(y) - 2)) * (1/length(x) + 1/length(y)))
        dft <- length(x) + length(y) - 2
        probt1 <- pt(t1, dft, lower.tail = T)
        probt2 <- pt(t2, dft, lower.tail = F)
        ifelse(probt1 < alpha & probt2 < alpha, decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval can be rejected", 
            decis <- "The null hypothesis that the difference between the means exceeds the equivalence interval cannot be rejected")
    }
    title <- "Schuirmann's Test of the Equivalence of Two Independent Groups"
    ei <- (c(equivint))
    names(ei) <- c("equivalence interval")
    tstats <- c(t1, t2)
    dfs <- c(dft, dft)
    pvals <- c(probt1, probt2)
    names(tstats) <- c("t1", "t2")
    names(dfs) <- c("dft1", "dft2")
    names(pvals) <- c("p_t1", "p_t2")
    list(title, ei = ei, tstats = tstats, dfs = dfs, pvals = pvals, decis = decis)
} 
